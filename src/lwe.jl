struct LweParams
    n :: Int32
    alpha_min :: Float64
    alpha_max :: Float64
end


struct LweKey
    params :: LweParams
    key :: Array{Int32, 1}

    function LweKey(rng::AbstractRNG, params::LweParams)
        new(params, rand_uniform_int32(rng, Int(params.n))) # TODO: remove Int()
    end

    # extractions Ring Lwe . Lwe
    function LweKey(params::LweParams, tlwe_key) # sans doute un param supplémentaire
        @assert isa(tlwe_key, TLweKey) # (can't do it in the signature because it's declared later)
        N = tlwe_key.params.N
        k = tlwe_key.params.k
        @assert params.n == k*N

        key = Array{Int32, 1}(params.n)
        key .= tlwe_key.key.coefs[:] # TODO: use an approprtiate method

        new(params, key)
    end
end


struct LweSampleArray
    a :: AbstractArray # the n coefs of the mask (times dims)
    b :: AbstractArray
    current_variances :: AbstractArray # average noise of the sample (times dims)

    LweSampleArray(params::LweParams, dims...) = new(
        Array{Torus32}(params.n, dims...),
        Array{Torus32}(dims...),
        Array{Float64}(dims...))
    LweSampleArray(a_arr, b_arr, c_arr) = new(a_arr, b_arr, c_arr)
end

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


# Sums an array over the axis, preserving the type
function sum_int32(arr::Array{Int32}, axis)
    # TODO: sum() of an Int32 array produces an Int64 array.
    # Putting in this hack here for the time being.
    s = squeeze(sum(arr, 1), 1)
    signed.(trunc.(UInt32, unsigned.(s) .& 0xffffffff))
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

    # !!! to preserve the RNG call order for testing
    for i in 1:size(result.b, 1)
        result.b[i] = rand_gaussian_torus32(rng, messages[i], alpha)
        result.a[:,i] .= rand_uniform_torus32(rng, n)
    end
    #result.b .= rand_gaussian_torus32(rng, Int32(0), alpha, size(messages)...) + messages
    #result.a .= rand_uniform_torus32(rng, n, size(messages)...)

    result.b .+= sum_int32(result.a .* key.key, 1)
    result.current_variances .= ones(Float64, size(messages)...) * alpha^2
end


# This function encrypts a message by using key and a given noise value
function lweSymEncryptWithExternalNoise(
        rng::AbstractRNG,
        result::LweSampleArray, messages::Array{Torus32}, noises::Array{Float64},
        alpha::Float64, key::LweKey)

    @assert size(result) == size(messages)
    @assert size(result) == size(noises)

    n = key.params.n

    result.b .= messages + dtot32.(noises)
    result.a .= rand_uniform_torus32(rng, n, size(messages)...)
    result.b .+= sum_int32(result.a .* key.key, 1)

    result.current_variances .= ones(Float64, size(messages)...) * alpha^2
end


# This function computes the phase of sample by using key : phi = b - a.s
function lwePhase(sample::LweSampleArray, key::LweKey)
    sample.b - sum_int32(sample.a .* key.key, 1)
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
function lweNoiselessTrivial(result::LweSampleArray, mu::Torus32, params::LweParams)
    lweNoiselessTrivial(result, ones(Torus32, size(result)...) * mu, params)
end
function lweNoiselessTrivial(result::LweSampleArray, mus::Array{Torus32}, params::LweParams)
    @assert size(mus) == size(result)
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
    result.current_variances .+= p^2 * sample.current_variances
end


# result = result - p.sample
function lweSubMulTo(result::LweSampleArray, p::Int32, sample::LweSampleArray, params::LweParams)
    result.a .-= p * sample.a
    result.b .-= p * sample.b
    result.current_variances .+= p^2 * sample.current_variances
end


struct LweKeySwitchKey
    n :: Int32 # length of the input key: s'
    t :: Int32 # decomposition length
    basebit :: Int32 # log_2(base)
    base :: Int32 # decomposition base: a power of 2
    out_params :: LweParams # params of the output key s
    # these don't seem to be used anywhere
    #ks0_raw :: Array{LweSample, 1} # tableau qui contient tout les Lwe samples de taille nlbase
    #ks1_raw :: Array{LweSample, 2} # de taille nl  pointe vers un tableau ks0_raw dont les cases sont espaceés de base positions
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
            n::Int32, t::Int32, basebit::Int32,
            in_key::LweKey, out_key::LweKey)

        out_params = out_key.params

        base = 1 << basebit
        ks = LweSampleArray(out_params, Int(base), Int(t), Int(n)) # TODO: get rid of Int()

        alpha = out_key.params.alpha_min

        # chose a random vector of gaussian noises
        noises = rand_gaussian_float(rng, alpha, base - 1, t, n)

        # recenter the noises
        noises -= mean(noises)

        # generate the ks

        # term h=0 as trivial encryption of 0 (it will not be used in the KeySwitching)
        lweNoiselessTrivial(view(ks, 1, 1:t, 1:n), zeros(Torus32, t, n), out_key.params)

        # mess::Torus32 = (in_key.key[i] * Int32(h - 1)) * Int32(1 << (32 - j * basebit))
        hs = Torus32(2):Torus32(base)
        js = Torus32(1):Torus32(t)
        messages = (
            reshape(in_key.key, 1, 1, Int(n))
            .* (reshape(hs, Int(base-1), 1, 1) - Torus32(1))
            .* (Torus32(1) .<< (32 - reshape(js, 1, Int(t), 1) * basebit))
            )

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
    outer_ind = collect(1:length(ind))

    mask = filter_func.(ind_flat)
    outer_ind_filtered = outer_ind[mask]
    ind_flat_filtered = ind_flat[mask]

    flat_coords = ind_dim * (outer_ind_filtered - 1) + ind_flat_filtered
    flat_result = src_flat[flat_coords]
    flat_result
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
        n::Int32, t::Int32, basebit::Int32)

    base::Int32 = 1 << basebit # base=2 in [CGGI16]
    prec_offset::Int32 = 1 << (32 - (1 + basebit * t)) # precision
    mask::Int32 = base - 1

    # ai is of size n
    js = reshape(1:t, Int(t), 1, 1)
    ai = reshape(ai, 1, Int(n), :)
    aijs = ((ai + prec_offset) .>> (32 - js * basebit)) .& mask + 1

    # TODO: batch over ciphertext bits too
    for i in 1:length(result)
        sub_ks = LweSampleArray(take_filtered(ks, aijs[:,:,i], j -> j != 1))
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
