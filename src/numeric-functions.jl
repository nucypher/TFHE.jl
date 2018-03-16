#=
numeric_functions.h
numeric_functions.cpp
=#


const Torus32 = Int32

const _two31 = Int64(1) << 31 # 2^31
const _two32 = Int64(1) << 32 # 2^32
const _two32_double = Float64(_two32)
const _two31_double = Float64(_two31)


generator = MersenneTwister(0)


#uniform_int_distribution<Torus32> uniformTorus32_distrib(INT32_MIN, INT32_MAX);
function uniformTorus32_distrib(generator)
    rand(generator, Torus32)
end


#uniform_int_distribution<int32_t> uniformInt_distrib(INT_MIN, INT_MAX);
function uniformInt_distrib(generator)
    rand(generator, Int32)
end


function tfhe_random_generator_setSeed(values)
    for value in values
        srand(generator, value)
    end
end


# Gaussian sample centered in message, with standard deviation sigma
function gaussian32(message::Torus32, sigma::Float64)
    # Attention: all the implementation will use the stdev instead of the gaussian fourier param
    err = randn(generator) * sigma
    message + dtot32(err)
end


# Used to approximate the phase to the nearest message possible in the message space
# The constant Msize will indicate on which message space we are working (how many messages possible)
#
# "travailler sur 63 bits au lieu de 64, car dans nos cas pratiques, c'est plus précis"
function approxPhase(phase::Torus32, Msize::Int32)
    interv::UInt64 = ((UInt64(1) << 63) / Msize) * 2 # width of each intervall
    half_interval::Uint64 = interv / 2 # begin of the first intervall
    phase64::Uint64 = (UInt64(phase) << 32) + half_interval
    # floor to the nearest multiples of interv
    phase64 -= phase64 % interv
    # rescale to torus32
    convert(Torus32, phase64 >> 32)
end


# Used to approximate the phase to the nearest message possible in the message space
# The constant Msize will indicate on which message space we are working (how many messages possible)
#
# "travailler sur 63 bits au lieu de 64, car dans nos cas pratiques, c'est plus précis"
modSwitchFromTorus32(phase::Torus32, Msize::Int) = modSwitchFromTorus32(phase, Int32(Msize))
function modSwitchFromTorus32(phase::Torus32, Msize::Int32)
    interv::UInt64 = ((UInt64(1) << 63)/Msize) * UInt64(2) # width of each intervall
    half_interval::UInt64 = interv / 2 # begin of the first intervall
    phase64::UInt64 = (UInt64(unsigned(phase)) << 32) + half_interval
    # floor to the nearest multiples of interv
    trunc(Int32, signed(div(phase64, interv)))
end

# Used to approximate the phase to the nearest message possible in the message space
# The constant Msize will indicate on which message space we are working (how many messages possible)
#
# "travailler sur 63 bits au lieu de 64, car dans nos cas pratiques, c'est plus précis"
modSwitchToTorus32(mu::Int, Msize::Int) = modSwitchToTorus32(Int32(mu), Int32(Msize))
function modSwitchToTorus32(mu::Int32, Msize::Int32)
    interv::UInt64 = ((UInt64(1) << 63) / Msize) * UInt64(2) # width of each intervall
    phase64::UInt64 = unsigned(mu) * interv
    # floor to the nearest multiples of interv
    Torus32(signed(UInt32(phase64 >> 32)))
end


# from double to Torus32
function dtot32(d::Float64)
    trunc(Int32, (d - trunc(d)) * 2^32)
end

# from Torus32 to double
function t32tod(x::Torus32)
    convert(Float64, x) / _two32_double
end

