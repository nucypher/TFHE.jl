module TFHE
using Random
using FFTW
using DarkIntegers

include("numeric-functions.jl")

include("polynomials.jl")

include("lwe.jl")

include("tlwe.jl")

include("tgsw.jl")

include("lwe-bootstrapping.jl")

include("keys.jl")
export tfhe_key_pair
export tfhe_parameters
export TFHEEncryptedBit
export TFHESecretKey
export TFHECloudKey
export tfhe_encrypt_bit
export tfhe_decrypt_bit

include("boot-gates.jl")
export tfhe_gate_NAND
export tfhe_gate_OR
export tfhe_gate_AND
export tfhe_gate_XOR
export tfhe_gate_XNOR
export tfhe_gate_NOT
export tfhe_gate_CONSTANT
export tfhe_gate_NOR
export tfhe_gate_ANDNY
export tfhe_gate_ANDYN
export tfhe_gate_ORNY
export tfhe_gate_ORYN
export tfhe_gate_MUX

end
