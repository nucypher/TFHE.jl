module TFHE

using Compat

if Base.thisminor(VERSION) > v"0.6"
    using FFTW
    using Random
    using LinearAlgebra
else
    const mul! = A_mul_B!
end


include("numeric-functions.jl")

include("polynomials.jl")

include("lwe.jl")

include("tlwe.jl")

include("tgsw.jl")

include("lwe-bootstrapping.jl")

include("keys.jl")
export tfhe_key_pair
export tfhe_parameters
export TFHESecretKey
export TFHECloudKey
export tfhe_encrypt
export tfhe_decrypt
export empty_ciphertext

include("boot-gates.jl")
export tfhe_gate_XNOR!
export tfhe_gate_CONSTANT!
export tfhe_gate_MUX!
export tfhe_gate_AND!
export tfhe_gate_NOT!

end
