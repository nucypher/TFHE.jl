struct LweParams
    n :: Int
    alpha_min :: Float64
    alpha_max :: Float64
end


struct LweKey
    params :: LweParams
    key :: Array{Int32, 1}

    function LweKey(rng::AbstractRNG, params::LweParams)
        new(params, rand_uniform_int32(rng, params.n))
    end

    # extractions Ring Lwe . Lwe
    function LweKey(params::LweParams, tlwe_key) # sans doute un param supplémentaire
        @assert isa(tlwe_key, TLweKey) # (can't do it in the signature because it's declared later)
        N = tlwe_key.params.N
        k = tlwe_key.params.k
        @assert params.n == k*N

        key = tlwe_key.key.coefs[:] # TODO: use an approprtiate method

        new(params, key)
    end
end


struct LweSampleArray{T, U, V}
    a :: T # the n coefs of the mask (times dims)
    b :: U
    current_variances :: V # average noise of the sample (times dims)
end

LweSampleArray(params::LweParams, dims...) = LweSampleArray(
    Array{Torus32}(params.n, dims...),
    Array{Torus32}(dims...),
    Array{Float64}(dims...))


Base.size(arr::LweSampleArray, args...) = size(arr.b, args...)
Base.view(arr::LweSampleArray, ranges...) = LweSampleArray(
    view(arr.a, 1:size(arr.a, 1), ranges...),
    view(arr.b, ranges...),
    view(arr.current_variances, ranges...))
Base.length(arr::LweSampleArray) = length(arr.b)
Base.reshape(arr::LweSampleArray, dims...) = LweSampleArray(
    reshape(arr.a, size(arr.a, 1), dims...),
    reshape(arr.b, dims...),
    reshape(arr.current_variances, dims...))
Base.getindex(arr::LweSampleArray, ranges...) = view(arr, ranges...)


function vec_mul_mat(b::Array{Int32, 1}, a::Union{Array{Int32}, SubArray{Int32}})
    s = squeeze(sum(a .* b, 1), 1)
    # TODO: sum() of an Int32 array produces an Int64 array.
    # Putting in this hack here for the time being.
    # This behavior may change in later Julia versions,
    # then the conversion will not be necessary.
    # (Despite the auto-promotion and conversion,
    # it is still faster than matrix multiplication.)
    to_int32.(s)
end


#=
 * This function encrypts message by using key, with stdev alpha
 * The Lwe sample for the result must be allocated and initialized
 * (this means that the parameters are already in the result)
=#
function lweSymEncrypt(
        rng::AbstractRNG, result::LweSampleArray,
        messages::Array{Torus32}, alpha::Float64, key::LweKey)

    @assert size(result) == size(messages)

    n = key.params.n

    result.b .= rand_gaussian_torus32(rng, Int32(0), alpha, size(messages)...) + messages
    result.a .= rand_uniform_torus32(rng, n, size(messages)...)
    result.b .+= vec_mul_mat(key.key, result.a)
    result.current_variances .= alpha^2
end


# This function encrypts a message by using key and a given noise value
function lweSymEncryptWithExternalNoise(
        rng::AbstractRNG,
        result::LweSampleArray, messages::Array{Torus32}, noises::Array{Float64},
        alpha::Float64, key::LweKey)

    @assert size(result) == size(messages)
    @assert size(result) == size(noises)

    n = key.params.n

    result.b .= messages .+ dtot32.(noises)
    result.a .= rand_uniform_torus32(rng, n, size(messages)...)
    result.b .+= vec_mul_mat(key.key, result.a)
    result.current_variances .= alpha^2
end


# This function computes the phase of sample by using key : phi = b - a.s
function lwePhase(sample::LweSampleArray, key::LweKey)
    sample.b - vec_mul_mat(key.key, sample.a)
end


# Arithmetic operations on Lwe samples


# result = sample
function lweCopy(result::LweSampleArray, sample::LweSampleArray, params::LweParams)
    result.a .= sample.a
    result.b .= sample.b
    result.current_variances .= sample.current_variances
end


# result = -sample
function lweNegate(result::LweSampleArray, sample::LweSampleArray, params::LweParams)
    result.a .= -sample.a
    result.b .= -sample.b
    result.current_variances .= sample.current_variances
end


# result = (0,mu)
function lweNoiselessTrivial(result::LweSampleArray, mus::Union{Array{Torus32}, Torus32}, params::LweParams)
    result.a .= 0
    result.b .= mus
    result.current_variances .= 0
end


# result = result + sample
function lweAddTo(result::LweSampleArray, sample::LweSampleArray, params::LweParams)
    result.a .+= sample.a
    result.b .+= sample.b
    result.current_variances .+= sample.current_variances
end


# result = result - sample
function lweSubTo(result::LweSampleArray, sample::LweSampleArray, params::LweParams)
    result.a .-= sample.a
    result.b .-= sample.b
    result.current_variances .+= sample.current_variances
end


function lweSubAll(result::LweSampleArray, samples::LweSampleArray, params::LweParams)
    for i in 1:length(samples)
        result.a .-= samples.a[:,i]
        result.b .-= samples.b[i]
        result.current_variances .+= samples.current_variances[i]
    end
end


# result = result + p.sample
function lweAddMulTo(result::LweSampleArray, p::Int32, sample::LweSampleArray, params::LweParams)
    result.a .+= p * sample.a
    result.b .+= p * sample.b
    result.current_variances .+= p^2 .* sample.current_variances
end


# result = result - p.sample
function lweSubMulTo(result::LweSampleArray, p::Int32, sample::LweSampleArray, params::LweParams)
    result.a .-= p * sample.a
    result.b .-= p * sample.b
    result.current_variances .+= p^2 .* sample.current_variances
end


struct LweKeySwitchKey
    n :: Int # length of the input key: s'
    t :: Int # decomposition length
    basebit :: Int # log_2(base)
    base :: Int # decomposition base: a power of 2
    out_params :: LweParams # params of the output key s
    ks :: LweSampleArray # the keyswitch elements: a n.l.base matrix
    # de taille n pointe vers ks1 un tableau dont les cases sont espaceés de ell positions

    #=
    Create the key switching key:
     * normalize the error in the beginning
     * chose a random vector of gaussian noises (same size as ks)
     * recenter the noises
     * generate the ks by creating noiseless encryprions and then add the noise
    =#
    function LweKeySwitchKey(
            rng::AbstractRNG,
            n::Int, t::Int, basebit::Int,
            in_key::LweKey, out_key::LweKey)

        out_params = out_key.params

        base = 1 << basebit
        ks = LweSampleArray(out_params, base, t, n)

        alpha = out_key.params.alpha_min

        # chose a random vector of gaussian noises
        noises = rand_gaussian_float(rng, alpha, base - 1, t, n)

        # recenter the noises
        noises .-= mean(noises)

        # generate the ks

        # term h=0 as trivial encryption of 0 (it will not be used in the KeySwitching)
        lweNoiselessTrivial(view(ks, 1, 1:t, 1:n), zeros(Torus32, t, n), out_key.params)

        # mess::Torus32 = (in_key.key[i] * Int32(h - 1)) * Int32(1 << (32 - j * basebit))
        hs = Torus32(2):Torus32(base)
        js = Torus32(1):Torus32(t)

        r_key = reshape(in_key.key, 1, 1, n)
        r_hs = reshape(hs, base - 1, 1, 1)
        r_js = reshape(js, 1, t, 1)

        messages = @. r_key * (r_hs - Torus32(1)) * (Torus32(1) << (32 - r_js * basebit))

        lweSymEncryptWithExternalNoise(rng, view(ks, 2:base, 1:t, 1:n), messages, noises, alpha, out_key)

        new(n, t, basebit, base, out_params, ks)
    end
end


function take_filtered(src, ind, filter_func)
    # src: array (N, dims...)
    # indices: array (dims...) with elements in the range 1:N
    # returns: an array (prod(dims)) where each element corresponds to
    #          an element `ind[i,j,...]` with value != 1,
    #          and contains `src[ind[i,j,...], i, j, ...]`
    src_flat = reshape(src, length(src))
    ind_flat = reshape(ind, length(ind))

    ind_dim = size(src, 1)
    outer_ind = 1:length(ind)

    mask = filter_func.(ind_flat)
    outer_ind_filtered = outer_ind[mask]
    ind_flat_filtered = ind_flat[mask]

    flat_coords = @. ind_dim * (outer_ind_filtered - 1) + ind_flat_filtered
    src_flat[flat_coords]
end


#=
 * translates the message of the result sample by -sum(a[i].s[i]) where s is the secret
 * embedded in ks.
 * @param result the LWE sample to translate by -sum(ai.si).
 * @param ks The (n x t x base) key switching key
 *        ks[i][j][k] encodes k.s[i]/base^(j+1)
 * @param params The common LWE parameters of ks and result
 * @param ai The input torus array
 * @param n The size of the input key
 * @param t The precision of the keyswitch (technically, 1/2.base^t)
 * @param basebit Log_2 of base
=#
function lweKeySwitchTranslate_fromArray(result::LweSampleArray,
        ks::LweSampleArray, params::LweParams,
        ai::Array{Torus32, 2},
        n::Int, t::Int, basebit::Int)

    base = 1 << basebit # base=2 in [CGGI16]
    prec_offset = 1 << (32 - (1 + basebit * t)) # precision
    mask = base - 1

    # ai is of size n
    js = reshape(1:t, t, 1, 1)
    ai = reshape(ai, 1, n, :)
    aijs = @. ((ai + prec_offset) >> (32 - js * basebit)) & mask + 1

    # TODO: batch over ciphertext bits too
    for i in 1:length(result)
        sub_ks = take_filtered(ks, view(aijs, :, :, i), j -> j != 1)
        lweSubAll(result[i], sub_ks, params)
    end
end


#sample=(a',b')
function lweKeySwitch(result::LweSampleArray, ks::LweKeySwitchKey, sample::LweSampleArray)

    params = ks.out_params
    n = ks.n
    basebit = ks.basebit
    t = ks.t

    lweNoiselessTrivial(result, sample.b, params)
    lweKeySwitchTranslate_fromArray(result, ks.ks, params, sample.a, n, t, basebit)
end
