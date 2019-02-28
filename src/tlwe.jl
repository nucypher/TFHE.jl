struct TLweParams
    polynomial_degree :: Int # a power of 2: degree of the polynomials
    mask_size :: Int # number of polynomials in the mask

    function TLweParams(polynomial_degree::Int, mask_size::Int)
        new(polynomial_degree, mask_size)
    end
end


struct TLweKey
    params :: TLweParams
    key :: Array{IntPolynomial, 1} # the key (i.e k binary polynomials)

    function TLweKey(rng::AbstractRNG, params::TLweParams)
        key = [
            int_polynomial(rand_uniform_bool(rng, params.polynomial_degree))
            for i in 1:params.mask_size]
        new(params, key)
    end
end


# extractions Ring Lwe . Lwe
function extract_lwe_key(tlwe_key::TLweKey)
    tlwe_params = tlwe_key.params

    key = vcat([poly.coeffs for poly in tlwe_key.key]...)

    LweKey(LweParams(tlwe_params.polynomial_degree * tlwe_params.mask_size), key)
end


mutable struct TLweSample
    params :: TLweParams
    a :: Array{TorusPolynomial, 1} # array of length mask_size+1: mask + right term
    current_variance :: Float64 # avg variance of the sample

    TLweSample(params::TLweParams, a::Array{TorusPolynomial, 1}, cv::Float64) = new(params, a, cv)
end


mutable struct TransformedTLweSample
    params :: TLweParams
    a :: Array{TransformedTorusPolynomial, 1} # array of length mask_size+1: mask + right term
    current_variance :: Float64 # avg variance of the sample

    TransformedTLweSample(
            params::TLweParams, a::Array{TransformedTorusPolynomial, 1}, cv::Float64) =
        new(params, a, cv)
end


function tlwe_extract_sample(x::TLweSample)
    a = vcat([reverse_polynomial(p).coeffs for p in x.a[1:end-1]]...)
    b = x.a[end].coeffs[1]
    LweSample(LweParams(length(a)), a, b, 0.) # TODO: calculate the current variance
end


# create an homogeneous tlwe sample
function tlwe_encrypt_zero(rng::AbstractRNG, alpha::Float64, key::TLweKey)
    params = key.params
    polynomial_degree = params.polynomial_degree
    mask_size = params.mask_size

    a_part = [torus_polynomial(rand_uniform_torus32(rng, polynomial_degree)) for i in 1:mask_size]
    a_last = (
        int_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, polynomial_degree))
        + sum(key.key .* a_part))
    TLweSample(params, vcat(a_part..., a_last), alpha^2)
end


# result = (0,mu)
function tlwe_noiseless_trivial(mu::TorusPolynomial, params::TLweParams)
    a_part = [zero_torus_polynomial(params.polynomial_degree) for i in 1:params.mask_size]
    a_last = deepcopy(mu)
    TLweSample(params, vcat(a_part..., a_last), 0.)
end


Base.:+(x::TLweSample, y::TLweSample) =
    TLweSample(x.params, x.a .+ y.a, x.current_variance + y.current_variance)


Base.:-(x::TLweSample, y::TLweSample) =
    TLweSample(x.params, x.a .- y.a, x.current_variance + y.current_variance)


shift_polynomial(x::TLweSample, shift::Integer) =
    TLweSample(x.params, shift_polynomial.(x.a, shift), x.current_variance)


forward_transform(x::TLweSample) =
    TransformedTLweSample(x.params, forward_transform.(x.a), x.current_variance)


inverse_transform(x::TransformedTLweSample) =
    TLweSample(x.params, inverse_transform.(x.a), x.current_variance)


# TODO: how to compute the variance correctly?
Base.:+(x::TransformedTLweSample, y::TransformedTLweSample) =
    TransformedTLweSample(x.params, x.a .+ y.a, x.current_variance + y.current_variance)


# TODO: how to compute the variance correctly?
Base.:*(x::TransformedTLweSample, y::TransformedTorusPolynomial) =
    TransformedTLweSample(x.params, x.a .* y, x.current_variance)
