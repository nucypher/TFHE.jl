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
            int_polynomial(rand_uniform_bool(rng, params.polynomial_degree))
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

    TLweSample(a::Array{TorusPolynomial, 1}, cv::Float64) = new(a, cv)
end


mutable struct TransformedTLweSample
    a :: Array{TransformedTorusPolynomial, 1} # array of length mask_size+1: mask + right term
    current_variance :: Float64 # avg variance of the sample

    function TransformedTLweSample(params::TLweParams)
        a = [
            TransformedTorusPolynomial(zeros(Complex{Float64}, params.polynomial_degree ÷ 2))
            for i in 1:(params.mask_size+1)]
        new(a, 0.)
    end

    TransformedTLweSample(a, cv) = new(a, cv)
end


function tlwe_extract_sample(x::TLweSample, rparams::TLweParams)
    a = vcat([reverse_polynomial(p).coeffs for p in x.a[1:end-1]]...)
    b = x.a[end].coeffs[1]
    LweSample(a, b, 0.) # TODO: calculate the current variance
end


# create an homogeneous tlwe sample
function tlwe_encrypt_zero(rng::AbstractRNG, alpha::Float64, key::TLweKey, params::TLweParams)
    polynomial_degree = params.polynomial_degree
    mask_size = params.mask_size

    a_part = [torusPolynomialUniform(rng, polynomial_degree) for i in 1:mask_size]
    a_last = (
        int_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, polynomial_degree))
        + sum(key.key .* a_part))
    TLweSample(vcat(a_part..., a_last), alpha^2)
end


# result = (0,mu)
function tlwe_noiseless_trivial(mu::TorusPolynomial, params::TLweParams)
    a_part = [zero_torus_polynomial(params.polynomial_degree) for i in 1:params.mask_size]
    a_last = deepcopy(mu)
    TLweSample(vcat(a_part..., a_last), 0.)
end


Base.:+(x::TLweSample, y::TLweSample) =
    TLweSample(x.a .+ y.a, x.current_variance + y.current_variance)


Base.:-(x::TLweSample, y::TLweSample) =
    TLweSample(x.a .- y.a, x.current_variance + y.current_variance)


DarkIntegers.shift_polynomial(s::TLweSample, shift::Integer) =
    TLweSample(shift_polynomial.(s.a, shift), s.current_variance)


forward_transform(source::TLweSample) =
    TransformedTLweSample(forward_transform.(source.a), source.current_variance)


inverse_transform(source::TransformedTLweSample) =
    TLweSample(inverse_transform.(source.a), source.current_variance)


Base.:+(x::TransformedTLweSample, y::TransformedTLweSample) =
    TransformedTLweSample(x.a .+ y.a, x.current_variance + y.current_variance) # TODO: how to compute the variance correctly?

Base.:*(x::TransformedTLweSample, y::TransformedTorusPolynomial) =
    TransformedTLweSample(x.a .* y, x.current_variance) # TODO: how to compute the variance correctly?
