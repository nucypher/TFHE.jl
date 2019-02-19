module TFHE
using Random
using FFTW
using DarkIntegers
using LinearAlgebra

include("numeric-functions.jl")

include("polynomials.jl")

include("lwe.jl")

include("tlwe.jl")

include("tgsw.jl")

include("keyswitch.jl")
include("bootstrap.jl")

include("keys.jl")
export tfhe_key_pair
export tfhe_parameters
export TFHEEncryptedBit
export TFHESecretKey
export TFHECloudKey
export tfhe_encrypt_bit
export tfhe_decrypt_bit

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
