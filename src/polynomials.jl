const Torus32 = Int32


# This structure represents an integer polynomial modulo X^N+1
mutable struct IntPolynomial
    N :: Int32
    coefs :: Array{Int32, 1}
    IntPolynomial(N::Int32) = new(N, Array{Int32}(undef, N))
end

function new_IntPolynomial_array(nbelts::Int32, N::Int32)
    [IntPolynomial(N) for i in 1:nbelts]
end


# This structure represents an torus polynomial modulo X^N+1
mutable struct TorusPolynomial
    N :: Int32
    coefsT :: Array{Torus32, 1}
    TorusPolynomial(N::Int32) = new(N, Array{Torus32}(undef, N))
end


function new_TorusPolynomial_array(nbelts::Int, N::Int32)
    [TorusPolynomial(N) for i in 1:nbelts]
end


# J: Used as a mock for the FFT processor; will be removed during cleanup
struct FFTProc
    N :: Int32
    Ns2 :: Int32
    _2N :: Int32
    omegaxminus1 :: Array{Complex{Float64}, 1}

    FFTProc(N::Int32) = new(N, N / 2, N * 2, [exp(im * x * pi / N) for x in 0:(N/2-1)])
end



#=
 * This structure is used for FFT operations, and is a representation
 * over C of a polynomial in R[X]/X^N+1
 * This type is meant to be specialized, and all implementations of the structure must be compatible
 * (reinterpret_cast) with this one. Namely, they should contain at most 2 pointers
=#
mutable struct LagrangeHalfCPolynomial
    coefsC :: Array{Complex{Float64}, 1}
    proc :: FFTProc

    function LagrangeHalfCPolynomial(N::Int32)
        #@assert N == 1024
        new(Array{Complex{Float64}}(undef, div(N, 2)), FFTProc(N))
    end
end


function IntPolynomial_ifft(result::LagrangeHalfCPolynomial, p::IntPolynomial)
    res = result.coefsC
    a = p.coefs
    N = p.N

    rev_in = Array{Float64}(undef, 2*N)
    for i in 1:N
        rev_in[i] = a[i]/2.
        rev_in[N + i] = -rev_in[i]
    end

    rev_out = rfft(rev_in)

    for i in 0:(div(N, 2)-1)
        res[i+1] = rev_out[2*i+1+1]
        @assert abs(rev_out[2*i+1]) < 1e-20
    end
end


function TorusPolynomial_ifft(result::LagrangeHalfCPolynomial, p::TorusPolynomial)
    res = result.coefsC
    a = p.coefsT
    N = p.N

    _2pm33::Float64 = 1. / Float64(Int64(1)<<33)


    rev_in = Array{Float64}(undef, 2*N)
    for i in 1:N
        rev_in[i] = a[i] * _2pm33
        rev_in[N + i] = -rev_in[i]
    end
    rev_out = rfft(rev_in)

    for i in 0:(div(N, 2) - 1)
        res[i+1] = rev_out[2*i+1+1]
        @assert abs(rev_out[2*i+1]) < 1e-20
    end
end


function TorusPolynomial_fft(result::TorusPolynomial, p::LagrangeHalfCPolynomial)
    res = result.coefsT
    a = p.coefsC
    N = result.N

    _2p32::Float64 = Float64(Int64(1) << 32)
    _1sN::Float64 = 1. / N

    fw_in = Array{Complex{Float64}}(undef, N+1)
    for i in 0:(div(N, 2))
        fw_in[2*i+1] = 0
    end
    for i in 0:(div(N, 2)-1)
        fw_in[2*i+1+1] = a[i+1] # conj
    end

    fw_out = irfft(fw_in, 2 * N) * (2 * N)

    for i in 0:(N-1)
        # pas besoin du fmod... Torus32(int64_t(fmod(rev_out[i]*_1sN,1.)*_2p32));
        res[i+1] = trunc(Torus32, round(Int64, fw_out[i+1] * _1sN * _2p32) << 32 >> 32)
    end

    for i in 0:(N-1)
        @assert abs(fw_out[N+i+1] + fw_out[i+1]) < 1e-20
    end
end


function torusPolynomialAddMulRFFT(result::TorusPolynomial, poly1::IntPolynomial, poly2::TorusPolynomial)
    N = poly1.N
    tmp = [LagrangeHalfCPolynomial(N) for i in 1:3]
    tmpr = TorusPolynomial(N)
    IntPolynomial_ifft(tmp[1], poly1)
    TorusPolynomial_ifft(tmp[2], poly2)
    LagrangeHalfCPolynomialMul(tmp[3], tmp[1], tmp[2])
    TorusPolynomial_fft(tmpr, tmp[3])
    torusPolynomialAddTo(result, tmpr)
end


torusPolynomialAddMulR = torusPolynomialAddMulRFFT


#MISC OPERATIONS

# sets to zero
function LagrangeHalfCPolynomialClear(reps::LagrangeHalfCPolynomial)
    for i in 1:reps.proc.Ns2
        reps.coefsC[i] = 0
    end
end


# termwise multiplication in Lagrange space */
function LagrangeHalfCPolynomialMul(
        result::LagrangeHalfCPolynomial,
        a::LagrangeHalfCPolynomial,
        b::LagrangeHalfCPolynomial)

    Ns2 = result.proc.Ns2
    aa = a.coefsC
    bb = b.coefsC
    rr = result.coefsC
    for i in 0:(Ns2-1)
        rr[i+1] = aa[i+1]*bb[i+1]
    end
end


# termwise multiplication and addTo in Lagrange space
function LagrangeHalfCPolynomialAddMul(
        accum::LagrangeHalfCPolynomial,
        a::LagrangeHalfCPolynomial,
        b::LagrangeHalfCPolynomial)

    Ns2 = accum.proc.Ns2
    aa = a.coefsC
    bb = b.coefsC
    rr = accum.coefsC
    for i in 0:(Ns2-1)
        rr[i+1] += aa[i+1]*bb[i+1]
    end
end


# Torus polynomial functions

# TorusPolynomial = 0
function torusPolynomialClear(result::TorusPolynomial)
    for i in 1:result.N
        result.coefsT[i] = 0
    end
end

# TorusPolynomial = random
function torusPolynomialUniform(rng::AbstractRNG, result::TorusPolynomial)
    for i in 1:result.N
        result.coefsT[i] = rand_uniform_torus32(rng)
    end
end

# TorusPolynomial = TorusPolynomial
function torusPolynomialCopy(result::TorusPolynomial, sample::TorusPolynomial)
    for i in 1:result.N
        result.coefsT[i] = sample.coefsT[i]
    end
end

# TorusPolynomial += TorusPolynomial
function torusPolynomialAddTo(result::TorusPolynomial, poly2::TorusPolynomial)
    for i in 1:result.N
        result.coefsT[i] += poly2.coefsT[i]
    end
end


# result = (X^ai-1) * source
function torusPolynomialMulByXaiMinusOne(
        result::TorusPolynomial, ai::Int32, source::TorusPolynomial)

    N = source.N
    out = result.coefsT
    in_ = source.coefsT

    @assert (ai >= 0 && ai < 2 * N)

    if ai < N
        for i in 0:(ai-1) # sur que i-a<0
            out[i + 1] = -in_[i - ai + N + 1] - in_[i + 1]
        end

        for i in ai:(N-1) # sur que N>i-a>=0
            out[i + 1] = in_[i - ai + 1] - in_[i + 1]
        end
    else
        aa = ai - N
        for i in 0:(aa-1) # sur que i-a<0
            out[i + 1] = in_[i - aa + N + 1] - in_[i + 1]
        end
        for i in aa:(N-1) # sur que N>i-a>=0
            out[i + 1] = -in_[i - aa + 1] - in_[i + 1]
        end
    end
end


# result= X^{a}*source
function torusPolynomialMulByXai(
        result::TorusPolynomial, a::Int32, source::TorusPolynomial)

    N = source.N
    out = result.coefsT
    in_ = source.coefsT

    @assert (a >= 0 && a < 2 * N)

    if a < N
        for i in 0:(a-1) # sur que i-a<0
            out[i + 1] = -in_[i - a + N + 1]
        end

        for i in a:(N-1) # sur que N>i-a>=0
            out[i + 1] = in_[i - a + 1]
        end
    else
        aa = a - N
        for i in 0:(aa-1) # sur que i-a<0
            out[i + 1] = in_[i - aa + N + 1]
        end
        for i in aa:(N-1) # sur que N>i-a>=0
            out[i + 1] = -in_[i - aa + 1]
        end
    end
end
