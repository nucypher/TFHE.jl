IntPolynomial = Polynomial{Int32, N} where N
TorusPolynomial = Polynomial{Torus32, N} where N


int_polynomial(coeffs) = Polynomial(coeffs, negacyclic_modulus)
torus_polynomial(coeffs) = Polynomial(coeffs, negacyclic_modulus)


zero_torus_polynomial(len) = torus_polynomial(zeros(Torus32, len))


struct TransformedTorusPolynomial
    coeffs :: Array{Complex{Float64}, 1}
end


Base.:+(x::TransformedTorusPolynomial, y::TransformedTorusPolynomial) =
    TransformedTorusPolynomial(x.coeffs .+ y.coeffs)


Base.:*(x::TransformedTorusPolynomial, y::TransformedTorusPolynomial) =
    TransformedTorusPolynomial(x.coeffs .* y.coeffs)


Base.broadcastable(x::TransformedTorusPolynomial) = Ref(x)


"""
For the given p(x), calculates p(1/x) with the applcation of the corresponding modulus
(x^N - 1) for cyclic polynomials, (x^N+1) for negacyclic ones.
"""
function reverse_polynomial(p::Polynomial)
    new_coeffs = collect(reverse(p.coeffs))
    mul_by_monomial(Polynomial(new_coeffs, p.modulus), length(new_coeffs) + 1)
end


#=
Using tangent FFT for polynomial convolution (by multiplying them in the transformed space).
Sacrificing some readability for a big speed improvement.
=#


struct ForwardTransformPlan

    plan :: Plan
    coeffs :: Array{Complex{Float64}, 1}
    buffer :: Array{Complex{Float64}, 1}

    function ForwardTransformPlan(len::Int)
        @assert len % 2 == 0
        idx = collect(0:len÷2-1)
        coeffs = exp.((-2im * pi / len / 2) .* idx)
        buffer = Array{Complex{Float64}, 1}(undef, len ÷ 2)
        plan = plan_fft(buffer)
        new(plan, coeffs, buffer)
    end
end


struct InverseTransformPlan

    plan :: Plan
    coeffs :: Array{Complex{Float64}, 1}
    complex_buffer :: Array{Complex{Float64}, 1}
    int_buffer :: Array{Torus32, 1}

    function InverseTransformPlan(len::Int)
        @assert len % 2 == 0
        idx = collect(0:len÷2-1)
        coeffs = exp.((-2im * pi / len / 2) .* idx)
        complex_buffer = Array{Complex{Float64}, 1}(undef, len ÷ 2)
        int_buffer = Array{Torus32, 1}(undef, len)
        plan = plan_ifft(complex_buffer)
        new(plan, coeffs, complex_buffer, int_buffer)
    end
end


_forward_transform_plans = Dict{Int, ForwardTransformPlan}()
_inverse_transform_plans = Dict{Int, InverseTransformPlan}()


function get_forward_transform_plan(len::Int)
    if !haskey(_forward_transform_plans, len)
        p = ForwardTransformPlan(len)
        _forward_transform_plans[len] = p
        p
    else
        _forward_transform_plans[len]
    end
end


function get_inverse_transform_plan(len::Int)
    if !haskey(_inverse_transform_plans, len)
        p = InverseTransformPlan(len)
        _inverse_transform_plans[len] = p
        p
    else
        _inverse_transform_plans[len]
    end
end


function forward_transform(p::Union{IntPolynomial, TorusPolynomial})
    c = p.coeffs
    N = length(c)
    p = get_forward_transform_plan(N)
    p.buffer .= (c[1:N÷2] .- im .* c[N÷2+1:end]) .* p.coeffs
    TransformedTorusPolynomial(p.plan * p.buffer)
end


to_int32(x::Int64) = signed(trunc(UInt32, unsigned(x) & 0xffffffff))
to_int32(x::Float64) = to_int32(round(Int64, x))


function inverse_transform(x::TransformedTorusPolynomial)

    c = x.coeffs
    len = length(c)
    N = length(c)*2
    p = get_inverse_transform_plan(N)

    mul!(p.complex_buffer, p.plan, x.coeffs)
    p.complex_buffer .= conj.(p.complex_buffer) .* p.coeffs
    p.int_buffer[1:len] .= to_int32.(real.(p.complex_buffer))
    p.int_buffer[len+1:end] .= to_int32.(imag.(p.complex_buffer))

    torus_polynomial(copy(p.int_buffer))
end


"""
Multiply integer and torus polynomial using convolution in the transformed space.

Warning: current hardcoded transformation method (FFT) has limited precision,
so the coefficients of the integer polynomial cannot be too large (up to 11 bits),
if the coefficients of the torus polynomial use all 32 bits.
"""
function transformed_mul(x::IntPolynomial, y::TorusPolynomial)
    inverse_transform(forward_transform(x) * forward_transform(y))
end
