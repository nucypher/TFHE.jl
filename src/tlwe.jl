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
    key :: IntPolynomialArray # the key (i.e k binary polynomials)

    function TLweKey(rng::AbstractRNG, params::TLweParams)
        N = params.N
        k = params.k
        key = IntPolynomialArray(Int(N), Int(k)) # TODO: remove Int()
        arr = flat_coefs(key)
        arr .= rand_uniform_int32(rng, Int(N), Int(k)) # TODO: remove Int()
        new(params, key)
    end
end


mutable struct TLweSample
    a :: TorusPolynomialArray # array of length k+1: mask + right term
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
        a = TorusPolynomialArray(Int(params.N), Int(k + 1)) # TODO: get rid of Int()
        #b = a + k;
        current_variance = 0;

        new(a, 0, k)
    end
end


function new_TLweSample_array(nbelts::Int, params::TLweParams)
    [TLweSample(params) for i in 1:nbelts]
end


mutable struct TLweSampleFFT
    a :: LagrangeHalfCPolynomialArray # array of length k+1: mask + right term
    #b :: LagrangeHalfCPolynomial # alias of a[k] to get the right term
    current_variance :: Float64 # avg variance of the sample
    k :: Int32 # required during the destructor call...

    function TLweSampleFFT(params::TLweParams)
        # a is a table of k+1 polynomials, b is an alias for &a[k]
        k = params.k
        a = LagrangeHalfCPolynomialArray(Int(params.N), Int(k + 1)) # TODO: get rid of Int()
        new(a, 0., k)
    end
end


function tLweExtractLweSampleIndex(
        result::LweSample, x::TLweSample, index::Int32, params::LweParams, rparams::TLweParams)

    N = rparams.N
    k = rparams.k
    @assert params.n == k*N

    # TODO: get rid of Int()
    # TODO: use an appropriate method to get coefsT
    a_view = reshape(result.a, Int(N), Int(k))
    a_view[1:(index+1),:] = x.a.coefsT[(index+1):-1:1, 1:k]
    a_view[(index+2):N,:] = -x.a.coefsT[N:-1:(index+2), 1:k]
    result.b = x.a.coefsT[index+1, k+1]
end


function tLweExtractLweSample(result::LweSample, x::TLweSample, params::LweParams, rparams::TLweParams)
    tLweExtractLweSampleIndex(result, x, Int32(0), params, rparams)
end


# create an homogeneous tlwe sample
function tLweSymEncryptZero(rng::AbstractRNG, result::TLweSample, alpha::Float64, key::TLweKey)
    N = key.params.N
    k = key.params.k

    # TODO: use an appropriate method

    result.a.coefsT[:,k+1] .= rand_gaussian_torus32(rng, Int32(0), alpha, Int(N)) # TODO: get rid of Int()

    a_part = view(result.a, 1:k)
    tp_uniform!(rng, a_part)

    for i in 1:k
        tp_add_mul!(
            view(result.a, (k+1):(k+1)),
            view(key.key, i:i),
            view(result.a, i:i))
    end

    result.current_variance = alpha * alpha
end


# Arithmetic operations on TLwe samples

# result = sample
function tLweCopy(result::TLweSample, sample::TLweSample, params::TLweParams)
    k = params.k
    N = params.N
    result.a.coefsT .= sample.a.coefsT # TODO: use an appropriate method?
    result.current_variance = sample.current_variance
end


# result = (0,mu)
function tLweNoiselessTrivial(result::TLweSample, mu::TorusPolynomialArray, params::TLweParams)
    k = params.k
    @assert size(mu) == ()
    tp_clear!(result.a)
    result.a.coefsT[:,result.k+1] .= mu.coefsT[:,1] # TODO: wrap in a function?
    result.current_variance = 0.
end


# result = result + sample
function tLweAddTo(result::TLweSample, sample::TLweSample, params::TLweParams)
    k = params.k
    tp_add_to!(result.a, sample.a)
    result.current_variance += sample.current_variance
end


# mult externe de X^ai-1 par bki
function tLweMulByXaiMinusOne(result::TLweSample, ai::Int32, bk::TLweSample, params::TLweParams)
    tp_mul_by_xai_minus_one!(result.a, ai, bk.a)
end


# Computes the inverse FFT of the coefficients of the TLWE sample
function tLweToFFTConvert(result::TLweSampleFFT, source::TLweSample, params::TLweParams)
    k = params.k
    tp_ifft!(result.a, source.a)
    result.current_variance = source.current_variance
end

# Computes the FFT of the coefficients of the TLWEfft sample
function tLweFromFFTConvert(result::TLweSample, source::TLweSampleFFT, params::TLweParams)
    k = params.k
    tp_fft!(result.a, source.a)
    result.current_variance = source.current_variance
end

# Arithmetic operations on TLwe samples

# result = (0,0)
function tLweFFTClear(result::TLweSampleFFT, params::TLweParams)
    k = params.k
    lp_clear!(result.a)
    result.current_variance = 0.
end

# result = result + p*sample
function tLweFFTAddMulRTo(
        result::TLweSampleFFT, p::LagrangeHalfCPolynomialArray,
        sample::TLweSampleFFT, params::TLweParams)

    k = params.k

    p = reshape(p, 1)
    lp_add_mul!(result.a, p, sample.a)
    # result.current_variance += sample.current_variance;
    # TODO: how to compute the variance correctly?
end
