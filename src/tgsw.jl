struct TGswParams
    l :: Int # decomp length
    Bgbit :: Int # log_2(Bg)
    Bg :: Int # decomposition base (must be a power of 2)
    halfBg :: Torus32 # Bg/2
    maskMod :: Torus32 # Bg-1
    tlwe_params :: TLweParams # Params of each row
    kpl :: Int # number of rows = (k+1)*l
    h :: Array{Torus32, 1} # powers of Bgbit
    offset :: Torus32 # offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))

    function TGswParams(l::Int, Bgbit::Int, tlwe_params::TLweParams)

        Bg = 1 << Bgbit
        halfBg = div(Bg, 2)

        h = @. Torus32(1) << (32 - (1:l) * Bgbit) # 1/(Bg^(i+1)) as a Torus32

        # offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))
        offset = to_int32(sum(1 .<< (32 .- (1:l) .* Bgbit)) * halfBg)

        new(
            l,
            Bgbit,
            Bg,
            halfBg,
            Bg - 1,
            tlwe_params,
            (tlwe_params.k + 1) * l,
            h,
            offset,
            )
    end
end


struct TGswKey
    params :: TGswParams # the parameters
    tlwe_params :: TLweParams # the tlwe params of each rows
    tlwe_key :: TLweKey

    function TGswKey(rng::AbstractRNG, params::TGswParams)
        tlwe_key = TLweKey(rng, params.tlwe_params)
        new(
            params,
            params.tlwe_params,
            tlwe_key
            )
    end
end


struct TGswSampleArray
    samples :: TLweSampleArray # TLweSample* all_sample; (k+1)l TLwe Sample
    k :: Int
    l :: Int

    function TGswSampleArray(params::TGswParams, dims...)
        k = params.tlwe_params.k
        l = params.l
        samples = TLweSampleArray(params.tlwe_params, l, k + 1, dims...)
        new(samples, k, l)
    end
end


struct TGswSampleFFTArray
    samples :: TLweSampleFFTArray # TLweSample* all_sample; (k+1)l TLwe Sample
    k :: Int
    l :: Int

    # constructor content
    function TGswSampleFFTArray(params::TGswParams, dims...)
        k = params.tlwe_params.k
        l = params.l
        samples = TLweSampleFFTArray(params.tlwe_params, l, k + 1, dims...)
        new(samples, k, l)
    end

    function TGswSampleFFTArray(samples, k, l)
        new(samples, k, l)
    end
end


Base.view(arr::TGswSampleFFTArray, ranges...) = TGswSampleFFTArray(
    view(arr.samples, 1:arr.l, 1:(arr.k + 1), ranges...), arr.k, arr.l)
Base.size(arr::TGswSampleFFTArray) = size(arr.samples)[3:end]


# Result += mu*H, mu integer
function tGswAddMuIntH(result::TGswSampleArray, messages::Array{Int32, 1}, params::TGswParams)
    k = params.tlwe_params.k
    l = params.l
    h = params.h

    # compute result += H

    # returns an underlying coefsT of TorusPolynomialArray, with the total size
    # (N, k + 1 [from TLweSample], l, k + 1 [from TGswSample], n)
    # messages: (n,)
    # h: (l,)
    # TODO: use an appropriate method
    # TODO: not sure if it's possible to fully vectorize it
    for bloc in 1:(k+1)
        result.samples.a.coefsT[1, bloc, :, bloc, :] .+= (
            reshape(messages, 1, length(messages))
            .* reshape(h, l, 1))
    end
end


# Result = tGsw(0)
function tGswEncryptZero(rng::AbstractRNG, result::TGswSampleArray, alpha::Float64, key::TGswKey)
    rlkey = key.tlwe_key
    tLweSymEncryptZero(rng, result.samples, alpha, rlkey)
end


# encrypts a constant message
function tGswSymEncryptInt(
        rng::AbstractRNG, result::TGswSampleArray, messages::Array{Int32, 1},
        alpha::Float64, key::TGswKey)
    tGswEncryptZero(rng, result, alpha, key)
    tGswAddMuIntH(result, messages, key.params)
end


function tGswTorus32PolynomialDecompH(sample::TorusPolynomialArray, params::TGswParams)

    N = params.tlwe_params.N
    l = params.l
    k = params.tlwe_params.k
    Bgbit = params.Bgbit

    maskMod = params.maskMod
    halfBg = params.halfBg
    offset = params.offset

    result = IntPolynomialArray(N, l, size(sample)...)

    decal(p) = 32 - p * Bgbit

    ps = reshape(1:l, 1, l, 1, 1)
    sample_coefs = reshape(sample.coefsT, N, 1, size(sample)...)

    # do the decomposition
    @. result.coefs = ((sample_coefs + offset) >> decal(ps)) & maskMod - halfBg

    result
end


# For all the kpl TLWE samples composing the TGSW sample
# It computes the inverse FFT of the coefficients of the TLWE sample
function tGswToFFTConvert(result::TGswSampleFFTArray, source::TGswSampleArray, params::TGswParams)
    tLweToFFTConvert(result.samples, source.samples, params.tlwe_params)
end


# External product (*): accum = gsw (*) accum
function tGswFFTExternMulToTLwe(accum::TLweSampleArray, gsw::TGswSampleFFTArray, params::TGswParams)
    tlwe_params = params.tlwe_params
    k = tlwe_params.k
    l = params.l
    kpl = params.kpl
    N = tlwe_params.N

    # TODO attention, improve these new/delete...
    tmpa = TLweSampleFFTArray(tlwe_params, size(accum)...)

    # shape: l, k + 1, size(accum)...
    deca = tGswTorus32PolynomialDecompH(accum.a, params)

    decaFFT = LagrangeHalfCPolynomialArray(N, l, k + 1, size(accum)...) # fft version
    ip_ifft!(decaFFT, deca)

    tLweFFTClear(tmpa, tlwe_params)

    for i in 1:(k+1)
        for j in 1:l
            tLweFFTAddMulRTo(
                tmpa, # N/2 x k+1 x message_len
                reshape(view(decaFFT, j, i, 1:size(accum, 1)), 1, size(accum, 1)), # N/2 x message_len
                view(gsw.samples, j, i), # N/2 x k+1
                tlwe_params)
        end
    end

    tLweFromFFTConvert(accum, tmpa, tlwe_params)
end
