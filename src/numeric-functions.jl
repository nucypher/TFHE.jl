const Torus32 = Int32


function rand_uniform_int32(rng::AbstractRNG, dims...)
    #=
    if length(dims) == 0
        rand(rng, Int32(0):Int32(1))
    else
        result = Array{Int32}(dims...)
        for i in 1:length(result)
            result[i] = rand(rng, Int32(0):Int32(1))
        end
        result
    end
    =#
    rand(rng, Int32(0):Int32(1), dims...)
end


function rand_uniform_torus32(rng::AbstractRNG, dims...)
    # TODO: if dims == () (it happens), the return value is not an array -> type instability
    #       also, there's probably instability for arrays of different dims too.
    #       Naturally, it applies for all other rand_ functions.
    # TODO: compatibility with the reference, remove later

    if length(dims) == 0
        rand(rng, Torus32)
    else
        result = Array{Torus32}(dims...)
        for i in 1:length(result)
            result[i] = rand(rng, Torus32)
        end
        result
    end

    #rand(rng, Torus32, dims...)
end


function rand_gaussian_float(rng::AbstractRNG, sigma::Float64, dims...)
    #=
    if length(dims) == 0
        randn(rng) * sigma
    else
        result = Array{Float64}(dims...)
        for i in 1:length(result)
            result[i] = rand(rng)
        end
        result * sigma
    end
    =#
    randn(rng, dims...) * sigma
end


# Gaussian sample centered in message, with standard deviation sigma
function rand_gaussian_torus32(rng::AbstractRNG, message::Torus32, sigma::Float64, dims...)
    # Attention: all the implementation will use the stdev instead of the gaussian fourier param

    if length(dims) == 0
        err = randn(rng) * sigma
    else
        err = Array{Float64}(dims...)
        for i in 1:length(err)
            err[i] = randn(rng)
        end
        err *= sigma
    end
    message + dtot32.(err)

    #err = randn(rng, dims...) * sigma
    #message + dtot32.(err)
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
function dtot32(d)
    trunc.(Int32, (d - trunc.(d)) * 2^32)
end


