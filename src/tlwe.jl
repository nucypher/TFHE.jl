#=
tlwe.h
tlwe_functions.h
lwe.cpp
tlwe-fft-operations.cpp
tlwe-functions.cpp
tlwe.cpp
=#

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

    TLweKey(params::TLweParams) = new(params, new_IntPolynomial_array(params.k, params.N))
end


mutable struct TLweSample
    a :: Array{TorusPolynomial, 1} # array of length k+1: mask + right term
    #b :: TorusPolynomial # alias of a[k] to get the right term
    current_variance :: Float64 # avg variance of the sample
    k :: Int32

    function TLweSample(params::TLweParams)
        # Small change here:
        # a is a table of k+1 polynomials, b is an alias for &a[k]
        # like that, we can access all the coefficients as before:
        #   &sample.a[0],...,&sample.a[k-1]  and &sample.b
        # or we can also do it in a single for loop
        #   &sample.a[0],...,&sample.a[k]
        k = params.k
        a = new_TorusPolynomial_array(k + 1, params.N)
        #b = a + k;
        current_variance = 0;

        new(a, 0, k)
    end
end


function new_TLweSample_array(nbelts::Int, params::TLweParams)
    [TLweSample(params) for i in 1:nbelts]
end


mutable struct TLweSampleFFT
    a :: Array{LagrangeHalfCPolynomial, 1} # array of length k+1: mask + right term
    #b :: LagrangeHalfCPolynomial # alias of a[k] to get the right term
    current_variance :: Float64 # avg variance of the sample
    k :: Int32 # required during the destructor call...

    function TLweSampleFFT(params::TLweParams)
        # a is a table of k+1 polynomials, b is an alias for &a[k]
        k = params.k
        a = [LagrangeHalfCPolynomial(params.N) for i in 1:(k+1)]
        new(a, 0., k)
    end
end


function tLweExtractLweSampleIndex(
        result::LweSample, x::TLweSample, index::Int32, params::LweParams, rparams::TLweParams)

    N = rparams.N
    k = rparams.k
    @assert params.n == k*N

    for i in 0:(k-1)
        for j in 0:index
            result.a[i*N+j+1] = x.a[i+1].coefsT[index-j+1]
        end
        for j in (index+1):(N-1)
            result.a[i*N+j+1] = -x.a[i+1].coefsT[N+index-j+1]
        end
        result.b = x.a[x.k+1].coefsT[index+1]
    end
end


function tLweExtractLweSample(result::LweSample, x::TLweSample, params::LweParams, rparams::TLweParams)
    tLweExtractLweSampleIndex(result, x, Int32(0), params, rparams)
end


# extractions Ring Lwe . Lwe
function tLweExtractKey(result::LweKey, key::TLweKey) # sans doute un param supplémentaire
    N = key.params.N
    k = key.params.k
    @assert result.params.n == k*N

    for i in 0:(k-1)
        for j in 0:(N-1)
            result.key[i*N+j+1] = key.key[i+1].coefs[j+1]
        end
    end
end


# TLwe
function tLweKeyGen(result::TLweKey)
    N = result.params.N
    k = result.params.k
    for i in 0:(k-1)
        for j in 0:(N-1)
            result.key[i+1].coefs[j+1] = rand(generator, 0:1)
        end
    end
end


# create an homogeneous tlwe sample
function tLweSymEncryptZero(result::TLweSample, alpha::Float64, key::TLweKey)
    N = key.params.N
    k = key.params.k

    for j in 0:(N-1)
        result.a[result.k+1].coefsT[j+1] = gaussian32(Int32(0), alpha)
    end

    for i in 0:(k-1)
        torusPolynomialUniform(result.a[i+1])
        torusPolynomialAddMulR(result.a[result.k+1], key.key[i+1], result.a[i+1])
    end

    result.current_variance = alpha * alpha
end


function tLweSymEncrypt(result::TLweSample, message::TorusPolynomial, alpha::Float64, key::TLweKey)
    N = key.params.N

    tLweSymEncryptZero(result, alpha, key)

    for j in 0:(N-1)
        result.a[result.k+1].coefsT[j+1] += message.coefsT[j+1]
    end
end


# encrypts a constant message
function tLweSymEncryptT(result::TLweSample, message::Torus32, alpha::Float64, key::TLweKey)
    tLweSymEncryptZero(result, alpha, key)
    result.a[result.k+1].coefsT[0+1] += message
end


# This function computes the phase of sample by using key : phi = b - a.s
function tLwePhase(phase::TorusPolynomial, sample::TLweSample, key::TLweKey)
    k = key.params.k
    torusPolynomialCopy(phase, sample.a[sample.k+1]) # phi = b
    for i in 0:(k-1)
        torusPolynomialSubMulR(phase, key.key[i+1], sample.a[i+1])
    end
end


# This function computes the approximation of the phase
# à revoir, surtout le Msize
function tLweApproxPhase(message::TorusPolynomial, phase::TorusPolynomial, Msize::Int32, N::Int32)
    for i in 0:(N-1)
        message.coefsT[i+1] = approxPhase(phase.coefsT[i+1], Msize)
    end
end


function tLweSymDecrypt(result::TorusPolynomial, sample::TLweSample, key::TLweKey, Msize::Int32)
    tLwePhase(result, sample, key)
    tLweApproxPhase(result, result, Msize, key.params.N)
end


function tLweSymDecryptT(sample::TLweSample, key::TLweKey, Msize::Int32)
    phase = TorusPolynomial(key.params.N)

    tLwePhase(phase, sample, key)
    approxPhase(phase.coefsT[0+1], Msize)
end


# Arithmetic operations on TLwe samples

# result = (0,0)
function tLweClear(result::TLweSample, params::TLweParams)
    k = params.k

    for i in 0:(k-1)
        torusPolynomialClear(result.a[i+1])
    end

    torusPolynomialClear(result.a[result.k+1])
    result.current_variance = 0.
end


# result = sample
function tLweCopy(result::TLweSample, sample::TLweSample, params::TLweParams)
    k = params.k
    N = params.N

    for i in 0:k
        for j in 0:(N-1)
            result.a[i+1].coefsT[j+1] = sample.a[i+1].coefsT[j+1]
        end
    end

    result.current_variance = sample.current_variance
end


# result = (0,mu)
function tLweNoiselessTrivial(result::TLweSample, mu::TorusPolynomial, params::TLweParams)
    k = params.k
    for i in 0:(k-1)
        torusPolynomialClear(result.a[i+1])
    end
    torusPolynomialCopy(result.a[result.k+1], mu)
    result.current_variance = 0.
end


# result = (0,mu) where mu is constant
function tLweNoiselessTrivialT(result::TLweSample, mu::Torus32, params::TLweParams)
    k = params.k

    for i in 0:(k-1)
        torusPolynomialClear(result.a[i+1])
    end
    torusPolynomialClear(result.a[result.k+1])
    result.a[result.k+1].coefsT[0+1] = mu
    result.current_variance = 0.
end

# result = result + sample
function tLweAddTo(result::TLweSample, sample::TLweSample, params::TLweParams)
    k = params.k

    for i in 0:(k-1)
        torusPolynomialAddTo(result.a[i+1], sample.a[i+1])
    end
    torusPolynomialAddTo(result.a[result.k+1], sample.a[sample.k+1])
    result.current_variance += sample.current_variance
end


# result = result - sample
function tLweSubTo(result::TLweSample, sample::TLweSample, params::TLweParams)
    k = params.k

    for i in 0:(k-1)
        torusPolynomialSubTo(result.a[i+1], sample.a[i+1])
    end
    torusPolynomialSubTo(result.a[result.k+1], sample.a[sample.k+1])
    result.current_variance += sample.current_variance
end


# result = result + p.sample
function tLweAddMulTo(result::TLweSample, p::Int32, sample::TLweSample, params::TLweParams)
    k = params.k

    for i in 0:(k-1)
        torusPolynomialAddMulZTo(result.a[i+1], p, sample.a[i+1])
    end
    torusPolynomialAddMulZTo(result.a[result.k+1], p, sample.a[sample.k+1])
    result.current_variance += (p * p) * sample.current_variance
end


# result = result - p.sample
function tLweSubMulTo(result::TLweSample, p::Int32, sample::TLweSample, params::TLweParams)
    k = params.k
    for i in 0:(k-1)
        torusPolynomialSubMulZTo(result.a[i+1], p, sample.a[i+1])
    end
    torusPolynomialSubMulZTo(result.a[result.k+1], p, sample.a[sample.k+1])
    result.current_variance += (p * p) * sample.current_variance
end


# result = result + p.sample
function
tLweAddMulRTo(result::TLweSample, p::IntPolynomial, sample::TLweSample, params::TLweParams)
    k = params.k
    for i in 0:k
        torusPolynomialAddMulR(result.a[i+1], p, sample.a[i+1])
    end
    result.current_variance += intPolynomialNormSq2(p) * sample.current_variance
end


# mult externe de X^ai-1 par bki
function tLweMulByXaiMinusOne(result::TLweSample, ai::Int32, bk::TLweSample, params::TLweParams)
    k = params.k
    for i in 0:k
        torusPolynomialMulByXaiMinusOne(result.a[i+1], ai, bk.a[i+1])
    end
end


# result += (0,x)
function tLweAddTTo(result::TLweSample, pos::Int32, x::Torus32, params::TLweParams)
    result.a[pos+1].coefsT[0+1] += x
end


# result += p*(0,x)
function tLweAddRTTo(result::TLweSample, pos::Int32, p::IntPolynomial, x::Torus32, params::TLweParams)
    N = params.N
    for i in 0:(N-1)
        result.a[pos+1].coefsT[i+1] += p.coefs[i+1] * x
    end
end


# Computes the inverse FFT of the coefficients of the TLWE sample
function tLweToFFTConvert(result::TLweSampleFFT, source::TLweSample, params::TLweParams)
    k = params.k

    for i in 0:k
        TorusPolynomial_ifft(result.a[i+1], source.a[i+1])
    end
    result.current_variance = source.current_variance
end

# Computes the FFT of the coefficients of the TLWEfft sample
function tLweFromFFTConvert(result::TLweSample, source::TLweSampleFFT, params::TLweParams)
    k = params.k

    for i in 0:k
        TorusPolynomial_fft(result.a[i+1], source.a[i+1])
    end
    result.current_variance = source.current_variance
end

# Arithmetic operations on TLwe samples

# result = (0,0)
function tLweFFTClear(result::TLweSampleFFT, params::TLweParams)
    k = params.k
    for i in 0:k
        LagrangeHalfCPolynomialClear(result.a[i+1])
    end
    result.current_variance = 0.
end

# result = result + p*sample
function tLweFFTAddMulRTo(
        result::TLweSampleFFT, p::LagrangeHalfCPolynomial, sample::TLweSampleFFT, params::TLweParams)

    k = params.k

    for i in 0:k
        LagrangeHalfCPolynomialAddMul(result.a[i+1], p, sample.a[i+1])
    end
    # result.current_variance += sample.current_variance;
    # TODO: how to compute the variance correctly?
end
