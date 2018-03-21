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


struct TLweSampleArray
    a :: TorusPolynomialArray # array of length k+1: mask + right term
    #b :: TorusPolynomial # alias of a[k] to get the right term
    current_variances :: AbstractArray # avg variance of the sample
    k :: Int32

    function TLweSampleArray(params::TLweParams, dims...)
        # Small change here:
        # a is a table of k+1 polynomials, b is an alias for &a[k]
        # like that, we can access all the coefficients as before:
        #   &sample.a[0],...,&sample.a[k-1]  and &sample.b
        # or we can also do it in a single for loop
        #   &sample.a[0],...,&sample.a[k]
        k = params.k
        a = TorusPolynomialArray(Int(params.N), Int(k + 1), dims...) # TODO: get rid of Int()
        current_variances = zeros(Float64, dims...)

        new(a, current_variances, k)
    end

    TLweSampleArray(a, cv, k) = new(a, cv, k)
end


mutable struct TLweSampleFFTArray
    a :: LagrangeHalfCPolynomialArray # array of length k+1: mask + right term
    #b :: LagrangeHalfCPolynomial # alias of a[k] to get the right term
    current_variances :: AbstractArray # avg variance of the sample
    k :: Int32 # required during the destructor call...

    function TLweSampleFFTArray(params::TLweParams, dims...)
        # a is a table of k+1 polynomials, b is an alias for &a[k]
        k = params.k
        a = LagrangeHalfCPolynomialArray(Int(params.N), Int(k + 1), dims...) # TODO: get rid of Int()
        current_variances = zeros(Float64, dims...)
        new(a, current_variances, k)
    end

    TLweSampleFFTArray(a, cv, k) = new(a, cv, k)
end


Base.size(arr::TLweSampleArray) = size(arr.current_variances)
Base.size(arr::TLweSampleFFTArray) = size(arr.current_variances)

Base.view(arr::TLweSampleArray, ranges...) = TLweSampleArray(
    view(arr.a, 1:(arr.k+1), ranges...),
    view(arr.current_variances, ranges...),
    arr.k)
Base.view(arr::TLweSampleFFTArray, ranges...) = TLweSampleFFTArray(
    view(arr.a, 1:(arr.k+1), ranges...),
    view(arr.current_variances, ranges...),
    arr.k)

Base.reshape(arr::TLweSampleFFTArray, dims...) = TLweSampleFFTArray(
    reshape(arr.a, Int(arr.k + 1), dims...),
    reshape(arr.current_variances, dims...),
    arr.k)


function tLweExtractLweSampleIndex(
        result::LweSampleArray, x::TLweSampleArray, index::Int32, params::LweParams, rparams::TLweParams)

    N = rparams.N
    k = rparams.k
    @assert params.n == k*N

    # TODO: get rid of Int()
    # TODO: use an appropriate method to get coefsT
    a_view = reshape(result.a, Int(N), Int(k), size(result)...)
    a_view[1:(index+1),:,:] .= x.a.coefsT[(index+1):-1:1, 1:k, :]
    a_view[(index+2):N,:,:] .= -x.a.coefsT[N:-1:(index+2), 1:k,:]

    result.b .= x.a.coefsT[index+1, k+1, :]
end


function tLweExtractLweSample(result::LweSampleArray, x::TLweSampleArray, params::LweParams, rparams::TLweParams)
    tLweExtractLweSampleIndex(result, x, Int32(0), params, rparams)
end


# create an homogeneous tlwe sample
function tLweSymEncryptZero(rng::AbstractRNG, rng2, result::TLweSampleArray, alpha::Float64, key::TLweKey)
    N = key.params.N
    k = key.params.k

    # TODO: use an appropriate method

    result.a.coefsT[:,k+1,:,:,:] .= rand_gaussian_torus32(rng, Int32(0), alpha, Int(N), size(result)...) # TODO: get rid of Int()

    a_part = view(result.a, 1:k,
        1:size(result.a.coefsT,3), 1:size(result.a.coefsT,4), 1:size(result.a.coefsT,5))
    tp_uniform!(rng2, a_part)

    for i in 1:k
        tp_add_mul!(
            view(result.a, (k+1):(k+1), 1:size(result.a.coefsT,3), 1:size(result.a.coefsT,4), 1:size(result.a.coefsT,5)),
            view(key.key, i:i),
            view(result.a, i:i, 1:size(result.a.coefsT,3), 1:size(result.a.coefsT,4), 1:size(result.a.coefsT,5)))
    end

    result.current_variances .= alpha^2
end


# Arithmetic operations on TLwe samples

# result = sample
function tLweCopy(result::TLweSampleArray, sample::TLweSampleArray, params::TLweParams)
    result.a.coefsT .= sample.a.coefsT # TODO: use an appropriate method?
    result.current_variances .= sample.current_variances
end


# result = (0,mu)
function tLweNoiselessTrivial(result::TLweSampleArray, mu::TorusPolynomialArray, params::TLweParams)
    k = params.k
    tp_clear!(result.a)
    result.a.coefsT[:,result.k+1,:] .= mu.coefsT[:,:] # TODO: wrap in a function?
    result.current_variances .= 0.
end


# result = result + sample
function tLweAddTo(result::TLweSampleArray, sample::TLweSampleArray, params::TLweParams)
    k = params.k
    tp_add_to!(result.a, sample.a)
    result.current_variances .+= sample.current_variances
end


# mult externe de X^ai-1 par bki
function tLweMulByXaiMinusOne(result::TLweSampleArray, ai::Array{Int32}, bk::TLweSampleArray, params::TLweParams)
    tp_mul_by_xai_minus_one!(result.a, ai, bk.a)
end


# Computes the inverse FFT of the coefficients of the TLWE sample
function tLweToFFTConvert(result::TLweSampleFFTArray, source::TLweSampleArray, params::TLweParams)
    tp_ifft!(result.a, source.a)
    result.current_variances .= source.current_variances
end

# Computes the FFT of the coefficients of the TLWEfft sample
function tLweFromFFTConvert(result::TLweSampleArray, source::TLweSampleFFTArray, params::TLweParams)
    tp_fft!(result.a, source.a)
    result.current_variances .= source.current_variances
end

# Arithmetic operations on TLwe samples

# result = (0,0)
function tLweFFTClear(result::TLweSampleFFTArray, params::TLweParams)
    lp_clear!(result.a)
    result.current_variances .= 0.
end

# result = result + p*sample
function tLweFFTAddMulRTo(
        result::TLweSampleFFTArray, p::LagrangeHalfCPolynomialArray,
        sample::TLweSampleFFTArray, params::TLweParams)
    lp_add_mul!(result.a, p, sample.a)
    # result.current_variance += sample.current_variance;
    # TODO: how to compute the variance correctly?
end
