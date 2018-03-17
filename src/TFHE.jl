module TFHE

#=
tfhe_core.h
tfhe_garbage_collector.h
tfhe_generic_streams.h
tfhe_generic_templates.h
tfhe_io.h

tfhe.h
tfhe_garbage_collector.cpp
tfhe_generic_streams.cpp
tfhe_io.cpp
=#

include("numeric-functions.jl")
export tfhe_random_generator_setSeed

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
export tfhe_encrypt_bit!
export tfhe_decrypt_bit!

include("boot-gates.jl")
export tfhe_gate_XNOR!
export tfhe_gate_CONSTANT!
export tfhe_gate_MUX!
export tfhe_gate_AND!
export tfhe_gate_NOT!

end
