struct TLweParams
    polynomial_degree :: Int # a power of 2: degree of the polynomials
    mask_size :: Int # number of polynomials in the mask
    min_noise :: Float64 # minimal noise s.t. the sample is secure
    max_noise :: Float64 # maximal noise s.t. we can decrypt
    extracted_lweparams :: LweParams # lwe params if one extracts

    function TLweParams(polynomial_degree::Int, mask_size::Int, min_noise::Float64, max_noise::Float64)
        new(
            polynomial_degree, mask_size, min_noise, max_noise,
            LweParams(polynomial_degree * mask_size, min_noise, max_noise))
    end
end


struct TLweKey
    params :: TLweParams # the parameters of the key
    key :: Array{IntPolynomial, 1} # the key (i.e k binary polynomials)

    function TLweKey(rng::AbstractRNG, params::TLweParams)
        key = [
            int_polynomial(rand_uniform_int32(rng, params.polynomial_degree))
            for i in 1:params.mask_size]
        new(params, key)
    end
end


mutable struct TLweSample
    a :: Array{TorusPolynomial, 1} # array of length mask_size+1: mask + right term
    current_variance :: Float64 # avg variance of the sample

    function TLweSample(params::TLweParams)
        a = [
            torus_polynomial(zeros(Torus32, params.polynomial_degree))
            for i in 1:(params.mask_size+1)]
        new(a, 0)
    end

    TLweSample(a, cv) = new(a, cv)
end


mutable struct TLweSampleFFT
    a :: Array{LagrangeHalfCPolynomial, 1} # array of length mask_size+1: mask + right term
    current_variance :: Float64 # avg variance of the sample

    function TLweSampleFFT(params::TLweParams)
        a = [
            LagrangeHalfCPolynomial(zeros(Complex{Float64}, params.polynomial_degree ÷ 2))
            for i in 1:(params.mask_size+1)]
        new(a, 0.)
    end

    TLweSampleFFT(a, cv) = new(a, cv)
end


function tLweExtractLweSampleIndex(
        x::TLweSample, index::Int32, rparams::TLweParams)

    params = rparams.extracted_lweparams

    result = LweSample(params)

    polynomial_degree = rparams.polynomial_degree
    mask_size = rparams.mask_size

    for i in 0:(mask_size-1)
        for j in 0:index
            result.a[i*polynomial_degree+j+1] = x.a[i+1].coeffs[index-j+1]
        end
        for j in (index+1):(polynomial_degree-1)
            result.a[i*polynomial_degree+j+1] = -x.a[i+1].coeffs[polynomial_degree+index-j+1]
        end
        result.b = x.a[mask_size+1].coeffs[index+1]
    end

    result
end


function tLweExtractLweSample(x::TLweSample, rparams::TLweParams)
    tLweExtractLweSampleIndex(x, Int32(0), rparams)
end


# create an homogeneous tlwe sample
function tLweSymEncryptZero(rng::AbstractRNG, alpha::Float64, key::TLweKey, params::TLweParams)
    polynomial_degree = params.polynomial_degree
    mask_size = params.mask_size

    result = TLweSample(params)

    for j in 0:(polynomial_degree-1)
        result.a[mask_size+1].coeffs[j+1] = rand_gaussian_torus32(rng, Int32(0), alpha)
    end

    for i in 0:(mask_size-1)
        result.a[i+1] = torusPolynomialUniform(rng, polynomial_degree)
        result.a[mask_size+1] += key.key[i+1] * result.a[i+1]
    end

    result.current_variance = alpha * alpha
    result
end


# Arithmetic operations on TLwe samples


# result = (0,mu)
function tLweNoiselessTrivial(mu::TorusPolynomial, params::TLweParams)
    result = TLweSample(params)
    for i in 0:(params.mask_size-1)
        torusPolynomialClear(result.a[i+1])
    end
    result.a[params.mask_size+1] = deepcopy(mu)
    result.current_variance = 0.
    result
end


# result = result + sample
function Base.:+(result::TLweSample, sample::TLweSample)
    TLweSample(result.a .+ sample.a, result.current_variance + sample.current_variance)
end


# mult externe de X^ai-1 par bki
function tLweMulByXaiMinusOne(ai::Integer, bk::TLweSample, params::TLweParams)
    result = TLweSample(params)
    for i in 0:params.mask_size
        result.a[i+1] = torusPolynomialMulByXaiMinusOne(ai, bk.a[i+1])
    end
    result
end


# Computes the inverse FFT of the coefficients of the TLWE sample
function tLweToFFTConvert(source::TLweSample, params::TLweParams)

    result = TLweSampleFFT(params)

    for i in 0:params.mask_size
        result.a[i+1] = forward_transform(source.a[i+1])
    end
    result.current_variance = source.current_variance

    result
end

# Computes the FFT of the coefficients of the TLWEfft sample
function tLweFromFFTConvert(source::TLweSampleFFT, params::TLweParams)

    result = TLweSample(params)

    for i in 0:params.mask_size
        result.a[i+1] = inverse_transform(source.a[i+1])
    end
    result.current_variance = source.current_variance

    result
end

# Arithmetic operations on TLwe samples

# result = (0,0)
function zero_tlwe_fft(params::TLweParams)
    result = TLweSampleFFT(params)
    for i in 0:params.mask_size
        LagrangeHalfCPolynomialClear(result.a[i+1])
    end
    result.current_variance = 0.
    result
end


Base.:+(x::TLweSampleFFT, y::TLweSampleFFT) =
    TLweSampleFFT(x.a .+ y.a, x.current_variance + y.current_variance) # TODO: how to compute the variance correctly?

Base.:*(x::TLweSampleFFT, y::LagrangeHalfCPolynomial) =
    TLweSampleFFT(x.a .* y, x.current_variance) # TODO: how to compute the variance correctly?
