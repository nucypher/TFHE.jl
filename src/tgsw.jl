struct TGswParams
    l :: Int32 # decomp length
    Bgbit :: Int32 # log_2(Bg)
    Bg :: Int32 # decomposition base (must be a power of 2)
    halfBg :: Int32 # Bg/2
    maskMod :: UInt32 # Bg-1
    tlwe_params :: TLweParams # Params of each row
    kpl :: Int32 # number of rows = (k+1)*l
    h :: Array{Torus32, 1} # powers of Bgbit
    offset :: UInt32 # offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))

    function TGswParams(l::Int32, Bgbit::Int32, tlwe_params::TLweParams)

        h = Array{Torus32}(undef, l)
        for i in 0:(l-1)
            kk::Int32 = 32 - (i + 1) * Bgbit
            h[i+1] = 1 << kk # 1/(Bg^(i+1)) as a Torus32
        end

        # offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))
        temp1::UInt32 = 0
        for i in 0:(l-1)
            temp0::UInt32 = 1 << (32 - (i + 1) * Bgbit)
            temp1 += temp0
        end
        Bg = 1 << Bgbit
        halfBg = Bg / 2
        offset = temp1 * halfBg

        new(
            l,
            Bgbit,
            Bg,
            halfBg,
            Bg - 1,
            tlwe_params,
            convert(Int32, (tlwe_params.k + 1) * l),
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


struct TGswSample
    samples :: Array{TLweSample, 2}

    function TGswSample(params::TGswParams)
        k = params.tlwe_params.k
        l = params.l
        samples = [TLweSample(params.tlwe_params) for i in 1:((k + 1) * l)]
        new(reshape(samples, Int64(l), k + 1))
    end

    TGswSample(samples) = new(samples)
end


struct TGswSampleFFT
    samples:: Array{TLweSampleFFT, 2}

    function TGswSampleFFT(params::TGswParams)
        k = params.tlwe_params.k
        l = params.l
        samples = [TLweSampleFFT(params.tlwe_params) for i in 1:((k + 1) * l)]
        new(reshape(samples, Int64(l), k + 1))
    end

    TGswSampleFFT(samples) = new(samples)
end


# Result += mu*H, mu integer
function tGswAddMuIntH(result::TGswSample, message::Int32, params::TGswParams)
    k = params.tlwe_params.k
    l = params.l
    h = params.h

    result = deepcopy(result)

    # compute result += H
    for bloc in 0:k
        for i in 0:(l-1)
            result.samples[i+1, bloc+1].a[bloc+1] += message * h[i+1]
        end
    end

    result
end


# Result = tGsw(0)
function tGswEncryptZero(rng::AbstractRNG, alpha::Float64, key::TGswKey, params::TGswParams)
    rlkey = key.tlwe_key
    kpl = key.params.kpl

    result = TGswSample(params)

    for p in 0:(kpl-1)
        result.samples[p+1] = tLweSymEncryptZero(rng, alpha, rlkey, params.tlwe_params)
    end

    result
end


# encrypts a constant message
function tGswSymEncryptInt(rng::AbstractRNG, message::Int32, alpha::Float64, key::TGswKey)
    result = tGswEncryptZero(rng, alpha, key, key.params)
    tGswAddMuIntH(result, message, key.params)
end


function tGswTorus32PolynomialDecompH(sample::TorusPolynomial, params::TGswParams)

    N = params.tlwe_params.N
    l = params.l
    Bgbit = params.Bgbit
    buf = sample.coeffs

    maskMod = params.maskMod
    halfBg = params.halfBg
    offset = params.offset

    [int_polynomial(
        (((buf .+ signed(offset)) .>> Int32((32 - (p + 1) * Bgbit))) .& Int32(maskMod)) .- halfBg)
        for p in 0:(l-1)
    ]
end


# For all the kpl TLWE samples composing the TGSW sample
# It computes the inverse FFT of the coefficients of the TLWE sample
function tGswToFFTConvert(source::TGswSample, params::TGswParams)
    kpl = params.kpl

    result = TGswSampleFFT(params)

    for p in 0:(kpl-1)
        result.samples[p+1] = tLweToFFTConvert(source.samples[p+1], params.tlwe_params)
    end

    result
end


# External product (*): accum = gsw (*) accum
function tGswFFTExternMulToTLwe(accum::TLweSample, gsw::TGswSampleFFT, params::TGswParams)
    tlwe_params = params.tlwe_params
    k = tlwe_params.k
    l = params.l
    kpl = params.kpl
    N = tlwe_params.N

    deca = vcat([tGswTorus32PolynomialDecompH(accum.a[i+1], params) for i in 0:k]...)
    decaFFT = IntPolynomial_ifft.(deca)

    tmpa = zero_tlwe_fft(tlwe_params)
    for p in 0:(kpl-1)
        tmpa += gsw.samples[p+1] * decaFFT[p+1]
    end

    tLweFromFFTConvert(tmpa, tlwe_params)
end
