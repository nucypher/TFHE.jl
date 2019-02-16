const Torus32 = Int32


IntPolynomial = Polynomial{Int32}
TorusPolynomial = Polynomial{Torus32}

int_polynomial(coeffs) = IntPolynomial(coeffs, true)
torus_polynomial(coeffs) = TorusPolynomial(coeffs, true)


struct LagrangeHalfCPolynomial
    coeffs :: Array{Complex{Float64}, 1}
end


Base.broadcastable(x::LagrangeHalfCPolynomial) = Ref(x)


struct ForwardTransformPlan

    plan
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

    plan
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
    LagrangeHalfCPolynomial(p.plan * p.buffer)
end


function to_int32(x::Int64)
    signed(trunc(UInt32, unsigned(x) & 0xffffffff))
end
function to_int32(x::Float64)
    to_int32(round(Int64, x))
end


function inverse_transform(x::LagrangeHalfCPolynomial)

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


#MISC OPERATIONS

# sets to zero
function LagrangeHalfCPolynomialClear(reps::LagrangeHalfCPolynomial)
    reps.coeffs .= 0
end


Base.:+(x::LagrangeHalfCPolynomial, y::LagrangeHalfCPolynomial) =
    LagrangeHalfCPolynomial(x.coeffs .+ y.coeffs)

Base.:*(x::LagrangeHalfCPolynomial, y::LagrangeHalfCPolynomial) =
    LagrangeHalfCPolynomial(x.coeffs .* y.coeffs)


# Torus polynomial functions

# TorusPolynomial = 0
function torusPolynomialClear(result::TorusPolynomial)
    result.coeffs .= 0
end

# TorusPolynomial = random
function torusPolynomialUniform(rng::AbstractRNG, N)
    torus_polynomial(rand_uniform_torus32(rng, N))
end


# result = (X^ai-1) * source
function torusPolynomialMulByXaiMinusOne(ai::Int32, source::TorusPolynomial)
    shift_polynomial(source, ai) - source
end


# result= X^{a}*source
function torusPolynomialMulByXai(a::Int32, source::TorusPolynomial)
    shift_polynomial(source, a)
end
