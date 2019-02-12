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

        key = Array{Int32}(undef, params.n)
        for i in 0:(k-1)
            for j in 0:(N-1)
                key[i*N+j+1] = tlwe_key.key[i+1].coefs[j+1]
            end
        end

        new(params, key)
    end
end


mutable struct LweSample
    a :: Array{Torus32, 1} # the n coefs of the mask
    b :: Torus32
    current_variance :: Float64 # average noise of the sample

    LweSample(params::LweParams) = new(Array{Torus32}(undef, params.n), 0, 0.)
end

const TFHEEncryptedBit = LweSample



#=
 * This function encrypts message by using key, with stdev alpha
 * The Lwe sample for the result must be allocated and initialized
 * (this means that the parameters are already in the result)
=#
function lweSymEncrypt(rng::AbstractRNG, result::LweSample, message::Torus32, alpha::Float64, key::LweKey)
    n = key.params.n

    result.b = rand_gaussian_torus32(rng, message, alpha)
    for i in 0:(n-1)
        result.a[i+1] = rand_uniform_torus32(rng)
        result.b += result.a[i+1] * key.key[i+1]
    end
    result.current_variance = alpha * alpha
end



# This function encrypts a message by using key and a given noise value
function lweSymEncryptWithExternalNoise(
        rng::AbstractRNG,
        result::LweSample, message::Torus32, noise::Float64, alpha::Float64, key::LweKey)

    n = key.params.n

    result.b = message + dtot32(noise)

    for i in 0:(n-1)
        result.a[i+1] = rand_uniform_torus32(rng)
        result.b += result.a[i+1] * key.key[i+1]
    end
    result.current_variance = alpha * alpha
end


# This function computes the phase of sample by using key : phi = b - a.s
function lwePhase(sample::LweSample, key::LweKey)
    n = key.params.n
    axs::Torus32 = 0
    a = sample.a
    k = key.key

    for i in 0:(n-1)
       axs += a[i+1] * k[i+1]
    end
    sample.b - axs
end


# Arithmetic operations on Lwe samples


# result = sample
function lweCopy(result::LweSample, sample::LweSample, params::LweParams)
    n = params.n

    for i in 0:(n-1)
        result.a[i+1] = sample.a[i+1]
    end
    result.b = sample.b
    result.current_variance = sample.current_variance
end


# result = -sample
function lweNegate(result::LweSample, sample::LweSample, params::LweParams)
    n = params.n

    for i in 0:(n-1)
        result.a[i+1] = -sample.a[i+1]
    end
    result.b = -sample.b
    result.current_variance = sample.current_variance
end


# result = (0,mu)
function lweNoiselessTrivial(result::LweSample, mu::Torus32, params::LweParams)
    n = params.n

    for i in 0:(n-1)
        result.a[i+1] = 0
    end
    result.b = mu
    result.current_variance = 0.
end

# result = result + sample
function lweAddTo(result::LweSample, sample::LweSample, params::LweParams)
    n = params.n

    for i in 0:(n-1)
        result.a[i+1] += sample.a[i+1]
    end
    result.b += sample.b
    result.current_variance += sample.current_variance
end


# result = result - sample
function lweSubTo(result::LweSample, sample::LweSample, params::LweParams)
    n = params.n
    sa = sample.a
    ra = result.a

    for i in 0:(n-1)
        ra[i+1] -= sa[i+1]
    end

    result.b -= sample.b
    result.current_variance += sample.current_variance
end

# result = result + p.sample
function lweAddMulTo(result::LweSample, p::Int32, sample::LweSample, params::LweParams)
    n = params.n

    for i in 0:(n-1)
        result.a[i+1] += p*sample.a[i+1]
    end
    result.b += p*sample.b
    result.current_variance += (p*p)*sample.current_variance
end

# result = result - p.sample
function lweSubMulTo(result::LweSample, p::Int32, sample::LweSample, params::LweParams)
    n = params.n

    for i in 0:(n-1)
        result.a[i+1] -= p*sample.a[i+1]
    end
    result.b -= p*sample.b
    result.current_variance += (p*p)*sample.current_variance
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
    ks :: Array{LweSample, 3} # the keyswitch elements: a n.l.base matrix
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

        base::Int32 = 1 << basebit
        ks0_raw = [LweSample(out_params) for i in 1:(n * t * base)]
        ks = reshape(ks0_raw, Int64(base), Int64(t), Int64(n))

        alpha::Float64 = out_key.params.alpha_min
        sizeks::Int32 = n * t * (base - 1)
        #const int32_t n_out = out_key.params.n;

        err::Float64 = 0

        # chose a random vector of gaussian noises
        noise = Array{Float64}(undef, sizeks)
        for i in 0:(sizeks-1)
            noise[i+1] = rand_gaussian_float(rng, alpha)
            err += noise[i+1]
        end
        # recenter the noises
        err = err / sizeks
        for i in 0:(sizeks-1)
            noise[i+1] -= err
        end

        # generate the ks
        index :: Int32 = 0
        for i in 0:(n-1)
            for j in 0:(t-1)

                # term h=0 as trivial encryption of 0 (it will not be used in the KeySwitching)
                lweNoiselessTrivial(ks[0+1,j+1,i+1], Torus32(0), out_key.params)
                for h in 1:(base-1) # pas le terme en 0
                    mess::Torus32 = (in_key.key[i+1] * Int32(h)) * Int32(1 << (32 - (j + 1) * basebit))
                    lweSymEncryptWithExternalNoise(
                        rng, ks[h+1,j+1,i+1], mess, noise[index+1], alpha, out_key)
                    index += 1
                end
            end
        end

        new(n, t, basebit, base, out_params, ks)
    end
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
function lweKeySwitchTranslate_fromArray(result::LweSample,
        ks::Array{LweSample, 3}, params::LweParams,
        ai::Array{Torus32, 1},
        n::Int32, t::Int32, basebit::Int32)

    base::Int32 = 1 << basebit # base=2 in [CGGI16]
    prec_offset::Int32 = 1 << (32 - (1 + basebit * t)) # precision
    mask::Int32 = base - 1

    for i in 0:(n-1)
        aibar::UInt32 = unsigned(ai[i+1] + prec_offset)
        for j in 0:(t-1)
            aij::UInt32 = (aibar >> (32 - (j + 1) * basebit)) & mask
            if aij != 0
                lweSubTo(result, ks[aij+1,j+1,i+1], params)
            end
        end
    end
end


#sample=(a',b')
function lweKeySwitch(result::LweSample, ks::LweKeySwitchKey, sample::LweSample)
    params = ks.out_params
    n = ks.n
    basebit = ks.basebit
    t = ks.t

    lweNoiselessTrivial(result, sample.b, params)
    lweKeySwitchTranslate_fromArray(result, ks.ks, params, sample.a, n, t, basebit)
end
