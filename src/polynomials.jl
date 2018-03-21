# This structure represents an integer polynomial modulo X^N+1
struct IntPolynomialArray
    coefs :: AbstractArray # Array{Int32}

    IntPolynomialArray(N::Int, dims...) = new(Array{Int32}(N, dims...))
    IntPolynomialArray(arr::AbstractArray) = new(arr)
end


# This structure represents an torus polynomial modulo X^N+1
struct TorusPolynomialArray
    coefsT :: AbstractArray # Array{Torus32}

    TorusPolynomialArray(N::Int, dims...) = new(Array{Torus32}(N, dims...))
    TorusPolynomialArray(arr::AbstractArray) = new(arr)
end


# This structure is used for FFT operations, and is a representation
# over C of a polynomial in R[X]/X^N+1
struct LagrangeHalfCPolynomialArray
    coefsC :: AbstractArray # Array{Complex{Float64}}

    function LagrangeHalfCPolynomialArray(N::Int, dims...)
        @assert mod(N, 2) == 0
        new(Array{Complex{Float64}}(div(N, 2), dims...))
    end

    LagrangeHalfCPolynomialArray(arr::AbstractArray) = new(arr)
end


flat_coefs(arr::IntPolynomialArray) = reshape(arr.coefs, size(arr.coefs, 1), length(arr))
flat_coefs(arr::TorusPolynomialArray) = reshape(arr.coefsT, size(arr.coefsT, 1), length(arr))
flat_coefs(arr::LagrangeHalfCPolynomialArray) = reshape(arr.coefsC, size(arr.coefsC, 1), length(arr))

polynomial_size(arr::IntPolynomialArray) = size(arr.coefs, 1)
polynomial_size(arr::TorusPolynomialArray) = size(arr.coefsT, 1)
polynomial_size(arr::LagrangeHalfCPolynomialArray) = size(arr.coefsC, 1) * 2

Base.size(arr::IntPolynomialArray) = size(arr.coefs)[2:end]
Base.size(arr::TorusPolynomialArray) = size(arr.coefsT)[2:end]
Base.size(arr::LagrangeHalfCPolynomialArray) = size(arr.coefsC)[2:end]

Base.length(arr::IntPolynomialArray) = prod(size(arr))
Base.length(arr::TorusPolynomialArray) = prod(size(arr))
Base.length(arr::LagrangeHalfCPolynomialArray) = prod(size(arr))

Base.reshape(arr::LagrangeHalfCPolynomialArray, dims...) =
    LagrangeHalfCPolynomialArray(reshape(arr.coefsC, size(arr.coefsC, 1), dims...))

Base.view(arr::IntPolynomialArray, ranges...) =
    IntPolynomialArray(view(arr.coefs, 1:size(arr.coefs, 1), ranges...))
Base.view(arr::TorusPolynomialArray, ranges...) =
    TorusPolynomialArray(view(arr.coefsT, 1:size(arr.coefsT, 1), ranges...))
Base.view(arr::LagrangeHalfCPolynomialArray, ranges...) =
    LagrangeHalfCPolynomialArray(view(arr.coefsC, 1:size(arr.coefsC, 1), ranges...))


function ip_ifft!(result::LagrangeHalfCPolynomialArray, p::IntPolynomialArray)
    res = flat_coefs(result)
    a = flat_coefs(p)
    N = polynomial_size(p)

    rev_in = Array{Float64, 2}(2 * N, length(p))
    rev_in[1:N,:] .= a / 2
    rev_in[N+1:end,:] .= -rev_in[1:N,:]

    # TODO: use a preallocated array and plan_rfft()
    rev_out = rfft(rev_in, 1)

    res[1:div(N,2),:] .= rev_out[2:2:N+1,:]
end


function tp_ifft!(result::LagrangeHalfCPolynomialArray, p::TorusPolynomialArray)
    res = flat_coefs(result)
    a = flat_coefs(p)
    N = polynomial_size(p)

    rev_in = Array{Float64, 2}(2 * N, length(p))
    rev_in[1:N,:] .= a[1:N,:] / 2^33
    rev_in[(N+1):end,:] .= -rev_in[1:N,:]

    # TODO: use a preallocated array and plan_rfft()
    rev_out = rfft(rev_in, 1)

    res[1:div(N,2),:] .= rev_out[2:2:N,:]
end


function tp_fft!(result::TorusPolynomialArray, p::LagrangeHalfCPolynomialArray)
    res = flat_coefs(result)
    a = flat_coefs(p)
    N = polynomial_size(p)

    fw_in = Array{Complex{Float64}, 2}(N + 1, length(p))
    fw_in[1:2:N+1,:] .= 0
    fw_in[2:2:N+1,:] .= a[:,:]

    # TODO: use a preallocated array and plan_irfft()
    fw_out = irfft(fw_in, 2 * N, 1) * (2 * N)

    # TODO: move to numeric-functions.jl (need to figure out how to preserve the broadcasting)
    # TODO: a view() is necessary here

    coeff = 2^32 / N
    res .= to_int32.(fw_out[1:N, :] .* coeff)
end


function tp_add_mul!(
        result::TorusPolynomialArray, poly1::IntPolynomialArray, poly2::TorusPolynomialArray)

    N = polynomial_size(result)
    tmp = [LagrangeHalfCPolynomialArray(N, size(result)...) for i in 1:3]
    tmpr = TorusPolynomialArray(N, size(result)...)

    ip_ifft!(tmp[1], poly1)
    tp_ifft!(tmp[2], poly2)
    lp_mul!(tmp[3], tmp[1], tmp[2])
    tp_fft!(tmpr, tmp[3])
    tp_add_to!(result, tmpr)
end


#MISC OPERATIONS

# sets to zero
function lp_clear!(reps::LagrangeHalfCPolynomialArray)
    reps.coefsC .= 0
end


# termwise multiplication in Lagrange space */
function lp_mul!(
        result::LagrangeHalfCPolynomialArray,
        a::LagrangeHalfCPolynomialArray,
        b::LagrangeHalfCPolynomialArray)

    result.coefsC .= a.coefsC .* b.coefsC
end


# termwise multiplication and addTo in Lagrange space
function lp_add_mul!(
        accum::LagrangeHalfCPolynomialArray,
        a::LagrangeHalfCPolynomialArray,
        b::LagrangeHalfCPolynomialArray)

    accum.coefsC .+= a.coefsC .* b.coefsC
end


# Torus polynomial functions

# TorusPolynomial = 0
function tp_clear!(result::TorusPolynomialArray)
    result.coefsT .= 0
end

# TorusPolynomial = random
function tp_uniform!(rng::AbstractRNG, result::TorusPolynomialArray)
    result.coefsT .= rand_uniform_torus32(rng, size(result.coefsT)...)
end

# TorusPolynomial += TorusPolynomial
function tp_add_to!(result::TorusPolynomialArray, poly2::TorusPolynomialArray)
    result.coefsT .+= poly2.coefsT
end


# result = (X^ai-1) * source
function tp_mul_by_xai_minus_one!(
        result::TorusPolynomialArray, ais::Array{Int32}, source::TorusPolynomialArray)

    N = polynomial_size(result)
    out = result.coefsT
    in_ = source.coefsT

    # TODO: we may need a fully vectorized form of this

    @assert (all(x >= 0 for x in ais) && all(x < 2 * N for x in ais))

    for i in 1:length(ais)
        ai = ais[i]
        if ai < N
            out[1:ai,:,i] .= -in_[(N-ai+1):N,:,i] - in_[1:ai,:,i] # sur que i-a<0
            out[(ai+1):N,:,i] .= in_[1:(N-ai),:,i] - in_[(ai+1):N,:,i] # sur que N>i-a>=0
        else
            aa = ai - N
            out[1:aa,:,i] .= in_[(N-aa+1):N,:,i] - in_[1:aa,:,i] # sur que i-a<0
            out[(aa+1):N,:,i] .= -in_[1:(N-aa),:,i] - in_[(aa+1):N,:,i] # sur que N>i-a>=0
        end
    end
end


# result= X^{a}*source
function tp_mul_by_xai!(
        result::TorusPolynomialArray, ais::Array{Int32}, source::TorusPolynomialArray)

    N = polynomial_size(result)
    out = result.coefsT
    in_ = source.coefsT

    # TODO: we may need a fully vectorized form of this

    @assert (all(x >= 0 for x in ais) && all(x < 2 * N for x in ais))

    for i in 1:length(ais)
        a = ais[i]
        if a < N
            out[1:a,i] .= -in_[(N-a+1):N,i] # sur que i-a<0
            out[(a+1):N,i] .= in_[1:(N-a),i] # sur que N>i-a>=0
        else
            aa = a - N
            out[1:aa,i] .= in_[(N-aa+1):N,i] # sur que i-a<0
            out[(aa+1):N,i] .= -in_[1:(N-aa),i] # sur que N>i-a>=0
        end
    end
end
