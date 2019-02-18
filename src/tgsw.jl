struct TGswParams
    decomp_length :: Int32 # decomposition length
    log2_base :: Int32 # log2(decomposition base)
    tlwe_params :: TLweParams # Params of each row

    base_powers :: Array{Torus32, 1} # powers of Bgbit
    offset :: UInt32 # offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))

    function TGswParams(l::Int, Bgbit::Int, tlwe_params::TLweParams)

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
            tlwe_params,
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
        mask_size = params.tlwe_params.mask_size
        l = params.decomp_length
        samples = [TLweSample(params.tlwe_params) for i in 1:((mask_size + 1) * l)]
        new(reshape(samples, Int64(l), mask_size + 1))
    end

    TGswSample(samples) = new(samples)
end


struct TGswSampleFFT
    samples:: Array{TransformedTLweSample, 2}

    function TGswSampleFFT(params::TGswParams)
        mask_size = params.tlwe_params.mask_size
        l = params.decomp_length
        samples = [TransformedTLweSample(params.tlwe_params) for i in 1:((mask_size + 1) * l)]
        new(reshape(samples, Int64(l), mask_size + 1))
    end

    TGswSampleFFT(samples) = new(samples)
end


# Result += mu*H, mu integer
function tGswAddMuIntH(result::TGswSample, message::Int32, params::TGswParams)
    mask_size = params.tlwe_params.mask_size
    l = params.decomp_length
    h = params.base_powers

    result = deepcopy(result)

    # compute result += H
    for bloc in 0:mask_size
        for i in 0:(l-1)
            result.samples[i+1, bloc+1].a[bloc+1] += message * h[i+1]
        end
    end

    result
end


# Result = tGsw(0)
function tGswEncryptZero(rng::AbstractRNG, alpha::Float64, key::TGswKey, params::TGswParams)
    rlkey = key.tlwe_key

    result = TGswSample(params)

    for p in 1:length(result.samples)
        result.samples[p] = tlwe_encrypt_zero(rng, alpha, rlkey, params.tlwe_params)
    end

    result
end


# encrypts a constant message
function tGswSymEncryptInt(rng::AbstractRNG, message::Int32, alpha::Float64, key::TGswKey)
    result = tGswEncryptZero(rng, alpha, key, key.params)
    tGswAddMuIntH(result, message, key.params)
end


function tGswTorus32PolynomialDecompH(sample::TorusPolynomial, params::TGswParams)

    l = params.decomp_length
    Bgbit = params.log2_base
    buf = sample.coeffs

    maskMod = 2^params.log2_base - 1
    halfBg = 2^(params.log2_base - 1)
    offset = params.offset

    [int_polynomial(
        (((buf .+ signed(offset)) .>> Int32((32 - (p + 1) * Bgbit))) .& Int32(maskMod)) .- halfBg)
        for p in 0:(l-1)
    ]
end


# For all the kpl TLWE samples composing the TGSW sample
# It computes the inverse FFT of the coefficients of the TLWE sample
function tGswToFFTConvert(source::TGswSample, params::TGswParams)

    result = TGswSampleFFT(params)

    for p in 1:length(result.samples)
        result.samples[p] = forward_transform(source.samples[p], params.tlwe_params)
    end

    result
end


# External product (*): accum = gsw (*) accum
function tGswFFTExternMulToTLwe(accum::TLweSample, gsw::TGswSampleFFT, params::TGswParams)
    tlwe_params = params.tlwe_params
    mask_size = tlwe_params.mask_size
    l = params.decomp_length

    deca = hcat([tGswTorus32PolynomialDecompH(accum.a[i+1], params) for i in 0:mask_size]...)
    decaFFT = forward_transform.(deca)

    tmpa = sum(gsw.samples .* decaFFT)

    inverse_transform(tmpa, tlwe_params)
end
