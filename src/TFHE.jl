module TFHE

using Random: AbstractRNG
using LinearAlgebra: mul!
using FFTW: plan_fft, plan_ifft, Plan
import DarkIntegers: shift_polynomial
using DarkIntegers: Polynomial

include("numeric-functions.jl")

include("polynomials.jl")

include("lwe.jl")

include("tlwe.jl")

include("tgsw.jl")

include("keyswitch.jl")

include("bootstrap.jl")

include("api.jl")
export make_key_pair
export LweSample
export SecretKey
export CloudKey
export encrypt
export decrypt

include("gates.jl")
export gate_nand
export gate_or
export gate_and
export gate_xor
export gate_xnor
export gate_not
export gate_constant
export gate_nor
export gate_andny
export gate_andyn
export gate_orny
export gate_oryn
export gate_mux

end
