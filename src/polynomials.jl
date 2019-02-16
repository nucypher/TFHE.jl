const Torus32 = Int32


IntPolynomial = Polynomial{Int32}
TorusPolynomial = Polynomial{Torus32}

int_polynomial(coeffs) = IntPolynomial(coeffs, true)
torus_polynomial(coeffs) = TorusPolynomial(coeffs, true)


struct LagrangeHalfCPolynomial
    coeffs :: Array{Complex{Float64}, 1}
end


Base.broadcastable(x::LagrangeHalfCPolynomial) = Ref(x)


function IntPolynomial_ifft(p::Union{IntPolynomial, TorusPolynomial})
    c = p.coeffs

    N = length(c)
    idx = collect(0:N÷2-1)
    fft_coeffs = exp.((-2im * pi / N / 2) .* idx)

    cc = c[1:N÷2] .- im .* c[N÷2+1:end]

    t = fft(cc .* fft_coeffs)

    LagrangeHalfCPolynomial(t)
end


function to_int32(x::Int64)
    signed(trunc(UInt32, unsigned(x) & 0xffffffff))
end
function to_int32(x::Float64)
    to_int32(round(Int64, x))
end


function TorusPolynomial_fft(p::LagrangeHalfCPolynomial)

    N = length(p.coeffs) * 2
    idx = collect(0:N÷2-1)
    fft_coeffs = exp.((-2im * pi / N / 2) .* idx)

    tt = conj.(ifft(p.coeffs)) .* fft_coeffs

    torus_polynomial(to_int32.(vcat(real.(tt), imag.(tt))))
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
