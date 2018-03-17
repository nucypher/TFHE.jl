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

        h = Array{Torus32, 1}(l)
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
    all_sample :: Array{TLweSample, 1} # TLweSample* all_sample; (k+1)l TLwe Sample
    bloc_sample :: Array{TLweSample, 2} # accès optionnel aux différents blocs de taille l.
    # double current_variance;
    k :: Int32
    l :: Int32

    # initialize the sample structure
    # (equivalent of the C++ constructor)
    function TGswSample(params::TGswParams)
        k = params.tlwe_params.k
        l = params.l
        # tous les samples comme un vecteur ligne
        all_sample = new_TLweSample_array((k + 1) * l, params.tlwe_params)
        # blocs horizontaux (l lignes) de la matrice TGsw
        bloc_sample = reshape(all_sample, Int64(l), k + 1)
        new(all_sample, bloc_sample, k, l)
    end

end


function new_TGswSample_array(nbelts::Int32, params::TGswParams)
    [TGswSample(params) for i in 1:nbelts]
end


struct TGswSampleFFT
    all_samples :: Array{TLweSampleFFT, 1} # TLweSample* all_sample; (k+1)l TLwe Sample
    sample:: Array{TLweSampleFFT, 2} # accès optionnel aux différents blocs de taille l.
    # double current_variance;
    k :: Int32
    l :: Int32

    # constructor content
    function TGswSampleFFT(params::TGswParams)
        k = params.tlwe_params.k
        l = params.l
        all_samples = [TLweSampleFFT(params.tlwe_params) for i in 1:((k + 1) * l)]
        sample = reshape(all_samples, Int64(l), k + 1)
        new(all_samples, sample, k, l)
    end
end


# Result += mu*H, mu integer
function tGswAddMuIntH(result::TGswSample, message::Int32, params::TGswParams)
    k = params.tlwe_params.k
    l = params.l
    h = params.h

    # compute result += H
    for bloc in 0:k
        for i in 0:(l-1)
            result.bloc_sample[i+1, bloc+1].a[bloc+1].coefsT[0+1] += message * h[i+1]
        end
    end
end


# Result = tGsw(0)
function tGswEncryptZero(rng::AbstractRNG, result::TGswSample, alpha::Float64, key::TGswKey)
    rlkey = key.tlwe_key
    kpl = key.params.kpl

    for p in 0:(kpl-1)
        tLweSymEncryptZero(rng, result.all_sample[p+1], alpha, rlkey)
    end
end


# encrypts a constant message
function tGswSymEncryptInt(rng::AbstractRNG, result::TGswSample, message::Int32, alpha::Float64, key::TGswKey)
    tGswEncryptZero(rng, result, alpha, key)
    tGswAddMuIntH(result, message, key.params)
end


function tGswTorus32PolynomialDecompH(sample::TorusPolynomial, params::TGswParams)

    N = params.tlwe_params.N
    l = params.l
    Bgbit = params.Bgbit
    buf = sample.coefsT

    maskMod = params.maskMod
    halfBg = params.halfBg
    offset = params.offset

    result = [IntPolynomial(N) for i in 1:l]

    # First, add offset to everyone
    for j in 0:(N-1)
        buf[j+1] += signed(offset)
    end

    # then, do the decomposition (in parallel)
    for p in 0:(l-1)
        decal::Int32 = (32 - (p + 1) * Bgbit)
        res_p = result[p+1].coefs
        for j in 0:(N-1)
            temp1::UInt32 = (buf[j+1] >> decal) & maskMod
            res_p[j+1] = temp1 - halfBg
        end
    end

    # finally, remove offset to everyone
    for j in 0:(N-1)
        buf[j+1] -= signed(offset)
    end

    result
end


# For all the kpl TLWE samples composing the TGSW sample
# It computes the inverse FFT of the coefficients of the TLWE sample
function tGswToFFTConvert(result::TGswSampleFFT, source::TGswSample, params::TGswParams)
    kpl = params.kpl

    for p in 0:(kpl-1)
        tLweToFFTConvert(result.all_samples[p+1], source.all_sample[p+1], params.tlwe_params)
    end
end


# External product (*): accum = gsw (*) accum
function tGswFFTExternMulToTLwe(accum::TLweSample, gsw::TGswSampleFFT, params::TGswParams)
    tlwe_params = params.tlwe_params
    k = tlwe_params.k
    l = params.l
    kpl = params.kpl
    N = tlwe_params.N

    # TODO attention, improve these new/delete...
    deca = new_IntPolynomial_array(kpl, N) # decomposed accumulator
    decaFFT = [LagrangeHalfCPolynomial(N) for i in 1:kpl] # fft version
    tmpa = TLweSampleFFT(tlwe_params)

    for i in 0:k
        deca[i * l + 1:(i+1)*l] = tGswTorus32PolynomialDecompH(accum.a[i+1], params)
    end

    for p in 0:(kpl-1)
        IntPolynomial_ifft(decaFFT[p+1], deca[p+1])
    end

    tLweFFTClear(tmpa, tlwe_params)

    for p in 0:(kpl-1)
        tLweFFTAddMulRTo(tmpa, decaFFT[p+1], gsw.all_samples[p+1], tlwe_params)
    end

    tLweFromFFTConvert(accum, tmpa, tlwe_params)
end
