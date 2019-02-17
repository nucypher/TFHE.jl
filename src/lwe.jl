struct LweParams
    len :: Int
    min_noise :: Float64
    max_noise :: Float64
end


struct LweKey
    params :: LweParams
    key :: Array{Int32, 1}

    function LweKey(rng::AbstractRNG, params::LweParams)
        new(params, rand_uniform_int32(rng, params.len))
    end

    # extractions Ring Lwe . Lwe
    function LweKey(params::LweParams, tlwe_key)
        @assert isa(tlwe_key, TLweKey) # (can't do it in the signature because it's declared later)
        polynomial_degree = tlwe_key.params.polynomial_degree
        mask_size = tlwe_key.params.mask_size
        @assert params.len == mask_size * polynomial_degree

        key = Array{Int32}(undef, params.len)
        for i in 0:(mask_size-1)
            for j in 0:(polynomial_degree-1)
                key[i*polynomial_degree+j+1] = tlwe_key.key[i+1].coeffs[j+1]
            end
        end

        new(params, key)
    end
end


mutable struct LweSample
    a :: Array{Torus32, 1} # the n coefs of the mask
    b :: Torus32
    current_variance :: Float64 # average noise of the sample

    LweSample(params::LweParams) = new(Array{Torus32}(undef, params.len), 0, 0.)
    LweSample(a, b, cv) = new(a, b, cv)
end

const TFHEEncryptedBit = LweSample



#=
 * This function encrypts message by using key, with stdev alpha
 * The Lwe sample for the result must be allocated and initialized
 * (this means that the parameters are already in the result)
=#
function lweSymEncrypt(rng::AbstractRNG, result::LweSample, message::Torus32, alpha::Float64, key::LweKey)
    n = key.params.len

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

    n = key.params.len

    result.b = message + dtot32(noise)

    for i in 0:(n-1)
        result.a[i+1] = rand_uniform_torus32(rng)
        result.b += result.a[i+1] * key.key[i+1]
    end
    result.current_variance = alpha * alpha
end


# This function computes the phase of sample by using key : phi = b - a.s
function lwePhase(sample::LweSample, key::LweKey)
    n = key.params.len
    axs::Torus32 = 0
    a = sample.a
    k = key.key

    for i in 0:(n-1)
       axs += a[i+1] * k[i+1]
    end
    sample.b - axs
end


# Arithmetic operations on Lwe samples


# result = -sample
function Base.:-(sample::LweSample)
    LweSample(-sample.a, -sample.b, sample.current_variance)
end


# result = (0,mu)
function lweNoiselessTrivial(mu::Torus32, params::LweParams)
    LweSample(zeros(Torus32, params.len), mu, 0.)
end


# result = result + sample
function Base.:+(result::LweSample, sample::LweSample)
    LweSample(
        result.a .+ sample.a, result.b + sample.b,
        result.current_variance + sample.current_variance)
end


# result = result - sample
function Base.:-(result::LweSample, sample::LweSample)
    LweSample(
        result.a .- sample.a, result.b - sample.b,
        result.current_variance + sample.current_variance)
end


function Base.:*(sample::LweSample, p::Integer)
    tp = Torus32(p)
    LweSample(sample.a .* tp, sample.b * tp, sample.current_variance * p^2)
end

Base.:*(p::Integer, sample::LweSample) = sample * p

