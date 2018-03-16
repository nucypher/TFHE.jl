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
export LweSample

include("tlwe.jl")

include("tgsw.jl")

include("lwe-bootstrapping.jl")

include("tfhe-gate-bootstrapping.jl")
export new_default_gate_bootstrapping_parameters
export new_random_gate_bootstrapping_secret_keyset
export new_gate_bootstrapping_ciphertext_array
export bootsSymEncrypt
export bootsSymDecrypt

export TFheGateBootstrappingCloudKeySet

include("boot-gates.jl")
export bootsXNOR
export bootsCONSTANT
export bootsMUX
export bootsAND
export bootsNOT


import Base: ==


function recursive_eq(x1, x2)
    if typeof(x1) != typeof(x2)
        return false
    end
    fs = fieldnames(x1)
    if length(fs) > 0
        all(recursive_eq(getfield(x1, f), getfield(x2, f)) for f in fs)
    else
        (x1 == x2)
    end
end

==(a::IntPolynomial, b::IntPolynomial) = recursive_eq(a, b)
==(a::TorusPolynomial, b::TorusPolynomial) = recursive_eq(a, b)
==(a::TGswSample, b::TGswSample) = recursive_eq(a, b)
==(a::TLweSample, b::TLweSample) = recursive_eq(a, b)
==(a::LweSample, b::LweSample) = recursive_eq(a, b)
==(a::LweBootstrappingKey, b::LweBootstrappingKey) = recursive_eq(a, b)
==(a::LweBootstrappingKeyFFT, b::LweBootstrappingKeyFFT) = recursive_eq(a, b)
==(a::TGswSampleFFT, b::TGswSampleFFT) = recursive_eq(a, b)
==(a::TLweSampleFFT, b::TLweSampleFFT) = recursive_eq(a, b)
==(a::LagrangeHalfCPolynomial, b::LagrangeHalfCPolynomial) = recursive_eq(a, b)


end
