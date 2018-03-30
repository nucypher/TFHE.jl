struct TLweParams
    N :: Int # a power of 2: degree of the polynomials
    k :: Int # number of polynomials in the mask
    alpha_min :: Float64 # minimal noise s.t. the sample is secure
    alpha_max :: Float64 # maximal noise s.t. we can decrypt
    extracted_lweparams :: LweParams # lwe params if one extracts

    function TLweParams(N::Int, k::Int, alpha_min::Float64, alpha_max::Float64)
        new(N, k, alpha_min, alpha_max, LweParams(N * k, alpha_min, alpha_max))
    end
end


struct TLweKey
    params :: TLweParams # the parameters of the key
    key :: IntPolynomialArray # the key (i.e k binary polynomials)

    function TLweKey(rng::AbstractRNG, params::TLweParams)
        N = params.N
        k = params.k
        key = IntPolynomialArray(rand_uniform_int32(rng, N, k))
        new(params, key)
    end
end


struct TLweSampleArray{T, U}
    a :: TorusPolynomialArray{T} # array of length k+1: mask + right term
    current_variances :: U # avg variance of the sample
    k :: Int
end


function TLweSampleArray(params::TLweParams, dims...)
    k = params.k
    a = TorusPolynomialArray(params.N, k + 1, dims...)
    current_variances = zeros(Float64, dims...)

    TLweSampleArray(a, current_variances, k)
end


mutable struct TLweSampleFFTArray{T, U}
    a :: LagrangeHalfCPolynomialArray{T} # array of length k+1: mask + right term
    current_variances :: U # avg variance of the sample
    k :: Int # required during the destructor call...
end


function TLweSampleFFTArray(params::TLweParams, dims...)
    # a is a table of k+1 polynomials
    k = params.k
    a = LagrangeHalfCPolynomialArray(params.N, k + 1, dims...)
    current_variances = zeros(Float64, dims...)
    TLweSampleFFTArray(a, current_variances, k)
end


Base.size(arr::TLweSampleArray, args...) = size(arr.current_variances, args...)
Base.size(arr::TLweSampleFFTArray, args...) = size(arr.current_variances, args...)


Base.reshape(arr::TLweSampleFFTArray, dims...) = TLweSampleFFTArray(
    reshape(arr.a, arr.k + 1, dims...),
    reshape(arr.current_variances, dims...),
    arr.k)


function tLweExtractLweSampleIndex(
        result::LweSampleArray, x::TLweSampleArray, index::Int, params::LweParams, rparams::TLweParams)

    N = rparams.N
    k = rparams.k
    @assert params.n == k*N

    # TODO: use an appropriate method to get coefsT
    a_view = reshape(result.a, N, k, size(result)...)
    a_view[1:(index+1),:,:] .= x.a.coefsT[(index+1):-1:1, 1:k, :]
    a_view[(index+2):N,:,:] .= .-x.a.coefsT[N:-1:(index+2), 1:k,:]

    result.b .= x.a.coefsT[index+1, k+1, :]

    nothing
end


function tLweExtractLweSample(result::LweSampleArray, x::TLweSampleArray, params::LweParams, rparams::TLweParams)
    tLweExtractLweSampleIndex(result, x, 0, params, rparams)
end


# create an homogeneous tlwe sample
function tLweSymEncryptZero(rng::AbstractRNG, result::TLweSampleArray, alpha::Float64, key::TLweKey)
    N = key.params.N
    k = key.params.k

    # TODO: use an appropriate method

    result.a.coefsT[:,k+1,:,:,:] .= rand_gaussian_torus32(rng, Int32(0), alpha, N, size(result)...)

    result.a.coefsT[:,1:k,:,:,:] .= rand_uniform_torus32(rng, N, k, size(result)...)

    tmp1 = LagrangeHalfCPolynomialArray(N, size(key.key)...)
    tmp2 = LagrangeHalfCPolynomialArray(N, k, size(result)...)
    tmp3 = LagrangeHalfCPolynomialArray(N, k, size(result)...)
    tmpr = TorusPolynomialArray(N, k, size(result)...)

    ip_ifft!(tmp1, key.key)
    tp_ifft!(tmp2, TorusPolynomialArray(result.a.coefsT[:, 1:k, :, :, :]))
    lp_mul!(tmp3, tmp1, tmp2)
    tp_fft!(tmpr, tmp3)

    for i in 1:k
        result.a.coefsT[:, k+1, :, :, :] .+= tmpr.coefsT[:, i, :, :, :]
    end

    result.current_variances .= alpha^2

    nothing
end


# Arithmetic operations on TLwe samples

# result = sample
function tLweCopy(result::TLweSampleArray, sample::TLweSampleArray, params::TLweParams)
    result.a.coefsT .= sample.a.coefsT # TODO: use an appropriate method?
    result.current_variances .= sample.current_variances
    nothing
end


# result = (0,mu)
function tLweNoiselessTrivial(result::TLweSampleArray, mu::TorusPolynomialArray, params::TLweParams)
    k = params.k
    tp_clear!(result.a)
    result.a.coefsT[:,result.k+1,:] .= mu.coefsT # TODO: wrap in a function?
    result.current_variances .= 0.
    nothing
end


# result = result + sample
function tLweAddTo(result::TLweSampleArray, sample::TLweSampleArray, params::TLweParams)
    k = params.k
    tp_add_to!(result.a, sample.a)
    result.current_variances .+= sample.current_variances
    nothing
end


# mult externe de X^ai-1 par bki
function tLweMulByXaiMinusOne(result::TLweSampleArray, ai::Array{Int32}, bk::TLweSampleArray, params::TLweParams)
    tp_mul_by_xai_minus_one!(result.a, ai, bk.a)
end


# Computes the inverse FFT of the coefficients of the TLWE sample
function tLweToFFTConvert(result::TLweSampleFFTArray, source::TLweSampleArray, params::TLweParams)
    tp_ifft!(result.a, source.a)
    result.current_variances .= source.current_variances
    nothing
end

# Computes the FFT of the coefficients of the TLWEfft sample
function tLweFromFFTConvert(result::TLweSampleArray, source::TLweSampleFFTArray, params::TLweParams)
    tp_fft!(result.a, source.a)
    result.current_variances .= source.current_variances
    nothing
end

# Arithmetic operations on TLwe samples

# result = (0,0)
function tLweFFTClear(result::TLweSampleFFTArray, params::TLweParams)
    lp_clear!(result.a)
    result.current_variances .= 0.
    nothing
end
