const Torus32 = Int32


function rand_uniform_int32(rng::AbstractRNG, dims...)
    rand(rng, Int32(0):Int32(1), dims...)
end


function rand_uniform_torus32(rng::AbstractRNG, dims...)
    # TODO: if dims == () (it happens), the return value is not an array -> type instability
    #       also, there's probably instability for arrays of different dims too.
    #       Naturally, it applies for all other rand_ functions.
    rand(rng, Torus32, dims...)
end


function rand_gaussian_float(rng::AbstractRNG, sigma::Float64, dims...)
    randn(rng, dims...) * sigma
end


# Gaussian sample centered in message, with standard deviation sigma
function rand_gaussian_torus32(rng::AbstractRNG, message::Torus32, sigma::Float64, dims...)
    # Attention: all the implementation will use the stdev instead of the gaussian fourier param
    err = randn(rng, dims...) * sigma
    message + dtot32.(err)
end


# Used to approximate the phase to the nearest message possible in the message space
# The constant Msize will indicate on which message space we are working (how many messages possible)
#
# "work on 63 bits instead of 64, because in our practical cases, it's more precise"
function modSwitchFromTorus32(phase::Torus32, Msize::Int)
    Msize = Int32(Msize)
    interv::UInt64 = ((UInt64(1) << 63)/Msize) * UInt64(2) # width of each intervall
    half_interval::UInt64 = interv / 2 # begin of the first intervall
    phase64::UInt64 = (UInt64(unsigned(phase)) << 32) + half_interval
    # floor to the nearest multiples of interv
    trunc(Int32, signed(div(phase64, interv)))
end

# Used to approximate the phase to the nearest message possible in the message space
# The constant Msize will indicate on which message space we are working (how many messages possible)
#
# "work on 63 bits instead of 64, because in our practical cases, it's more precise"
function modSwitchToTorus32(mu::Int, Msize::Int)
    mu = Int32(mu)
    Msize = Int32(Msize)
    interv::UInt64 = ((UInt64(1) << 63) / Msize) * UInt64(2) # width of each intervall
    phase64::UInt64 = unsigned(mu) * interv
    # floor to the nearest multiples of interv
    Torus32(signed(UInt32(phase64 >> 32)))
end


# from double to Torus32
function dtot32(d)
    trunc.(Int32, (d - trunc.(d)) * 2^32)
end


