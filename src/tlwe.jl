struct TLweParams
    N :: Int32 # a power of 2: degree of the polynomials
    k :: Int32 # number of polynomials in the mask
    alpha_min :: Float64 # minimal noise s.t. the sample is secure
    alpha_max :: Float64 # maximal noise s.t. we can decrypt
    extracted_lweparams :: LweParams # lwe params if one extracts

    function TLweParams(N::Int32, k::Int32, alpha_min::Float64, alpha_max::Float64)
        new(N, k, alpha_min, alpha_max, LweParams(N * k, alpha_min, alpha_max))
    end
end


struct TLweKey
    params :: TLweParams # the parameters of the key
    key :: Array{IntPolynomial, 1} # the key (i.e k binary polynomials)

    function TLweKey(rng::AbstractRNG, params::TLweParams)
        N = params.N
        k = params.k
        key = [int_polynomial(rand_uniform_int32(rng, N)) for i in 0:(k-1)]
        new(params, key)
    end
end


mutable struct TLweSample
    a :: Array{TorusPolynomial, 1} # array of length k+1: mask + right term
    #b :: TorusPolynomial # alias of a[k] to get the right term
    current_variance :: Float64 # avg variance of the sample

    function TLweSample(params::TLweParams)
        # Small change here:
        # a is a table of k+1 polynomials, b is an alias for &a[k]
        # like that, we can access all the coefficients as before:
        #   &sample.a[0],...,&sample.a[k-1]  and &sample.b
        # or we can also do it in a single for loop
        #   &sample.a[0],...,&sample.a[k]
        k = params.k
        a = [torus_polynomial(zeros(Torus32, params.N)) for i in 1:(k+1)]
        #b = a + k;
        current_variance = 0

        new(a, 0)
    end

    TLweSample(a, cv) = new(a, cv)
end


mutable struct TLweSampleFFT
    a :: Array{LagrangeHalfCPolynomial, 1} # array of length k+1: mask + right term
    #b :: LagrangeHalfCPolynomial # alias of a[k] to get the right term
    current_variance :: Float64 # avg variance of the sample

    function TLweSampleFFT(params::TLweParams)
        # a is a table of k+1 polynomials, b is an alias for &a[k]
        k = params.k
        a = [LagrangeHalfCPolynomial(zeros(Complex{Float64}, params.N ÷ 2)) for i in 1:(k+1)]
        new(a, 0.)
    end

    TLweSampleFFT(a, cv) = new(a, cv)
end


function tLweExtractLweSampleIndex(
        x::TLweSample, index::Int32, rparams::TLweParams)

    params = rparams.extracted_lweparams

    result = LweSample(params)

    N = rparams.N
    k = rparams.k
    @assert params.n == k*N

    for i in 0:(k-1)
        for j in 0:index
            result.a[i*N+j+1] = x.a[i+1].coeffs[index-j+1]
        end
        for j in (index+1):(N-1)
            result.a[i*N+j+1] = -x.a[i+1].coeffs[N+index-j+1]
        end
        result.b = x.a[k+1].coeffs[index+1]
    end

    result
end


function tLweExtractLweSample(x::TLweSample, rparams::TLweParams)
    tLweExtractLweSampleIndex(x, Int32(0), rparams)
end


# create an homogeneous tlwe sample
function tLweSymEncryptZero(rng::AbstractRNG, alpha::Float64, key::TLweKey, params::TLweParams)
    N = key.params.N
    k = key.params.k

    result = TLweSample(params)

    for j in 0:(N-1)
        result.a[k+1].coeffs[j+1] = rand_gaussian_torus32(rng, Int32(0), alpha)
    end

    for i in 0:(k-1)
        result.a[i+1] = torusPolynomialUniform(rng, N)
        result.a[k+1] += key.key[i+1] * result.a[i+1]
    end

    result.current_variance = alpha * alpha
    result
end


# Arithmetic operations on TLwe samples


# result = (0,mu)
function tLweNoiselessTrivial(mu::TorusPolynomial, params::TLweParams)
    result = TLweSample(params)
    k = params.k
    for i in 0:(k-1)
        torusPolynomialClear(result.a[i+1])
    end
    result.a[k+1] = deepcopy(mu)
    result.current_variance = 0.
    result
end


# result = result + sample
function Base.:+(result::TLweSample, sample::TLweSample)
    TLweSample(result.a .+ sample.a, result.current_variance + sample.current_variance)
end


# mult externe de X^ai-1 par bki
function tLweMulByXaiMinusOne(ai::Int32, bk::TLweSample, params::TLweParams)
    result = TLweSample(params)
    k = params.k
    for i in 0:k
        result.a[i+1] = torusPolynomialMulByXaiMinusOne(ai, bk.a[i+1])
    end
    result
end


# Computes the inverse FFT of the coefficients of the TLWE sample
function tLweToFFTConvert(source::TLweSample, params::TLweParams)

    result = TLweSampleFFT(params)

    k = params.k

    for i in 0:k
        result.a[i+1] = IntPolynomial_ifft(source.a[i+1])
    end
    result.current_variance = source.current_variance

    result
end

# Computes the FFT of the coefficients of the TLWEfft sample
function tLweFromFFTConvert(source::TLweSampleFFT, params::TLweParams)

    result = TLweSample(params)

    k = params.k

    for i in 0:k
        result.a[i+1] = TorusPolynomial_fft(source.a[i+1])
    end
    result.current_variance = source.current_variance

    result
end

# Arithmetic operations on TLwe samples

# result = (0,0)
function zero_tlwe_fft(params::TLweParams)
    result = TLweSampleFFT(params)
    k = params.k
    for i in 0:k
        LagrangeHalfCPolynomialClear(result.a[i+1])
    end
    result.current_variance = 0.
    result
end


Base.:+(x::TLweSampleFFT, y::TLweSampleFFT) =
    TLweSampleFFT(x.a .+ y.a, x.current_variance + y.current_variance) # TODO: how to compute the variance correctly?

Base.:*(x::TLweSampleFFT, y::LagrangeHalfCPolynomial) =
    TLweSampleFFT(x.a .* y, x.current_variance) # TODO: how to compute the variance correctly?
