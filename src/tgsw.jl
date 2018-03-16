#=
tgsw.h
tgsw_functions.h
tgsw-fft-operations.cpp
tgsw-functions.cpp
tgsw.cpp
=#

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
    key :: Array{IntPolynomial, 1} # the key (array of k polynomials)
    tlwe_key :: TLweKey

    function TGswKey(params::TGswParams)
        tlwe_key = TLweKey(params.tlwe_params)
        new(
            params,
            params.tlwe_params,
            tlwe_key.key,
            tlwe_key
            )
    end
end


# TGsw
# generate a tgsw key (in fact, a tlwe key)
function tGswKeyGen(result::TGswKey)
    tLweKeyGen(result.tlwe_key)
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


function tGswClear(result::TGswSample, params::TGswParams)
    kpl = params.kpl
    for p in 0:(kpl-1)
        tLweClear(result.all_sample[p+1], params.tlwe_params)
    end
end


# Result += H
function tGswAddH(result::TGswSample, params::TGswParams)
    k = params.tlwe_params.k
    l = params.l
    h = params.h

    # compute result += H
    for bloc in 0:k
        for i in 0:(l-1)
            result.bloc_sample[i+1, bloc+1].a[bloc+1].coefsT[0+1] += h[i+1]
        end
    end
end


# Result += mu*H
function tGswAddMuH(result::TGswSample, message::IntPolynomial, params::TGswParams)
    k = params.tlwe_params.k
    N = params.tlwe_params.N
    l = params.l
    h = params.h
    mu = message.coefs

    # compute result += H
    for bloc in 0:k
        for i in 0:(l-1)
            target = result.bloc_sample[i+1, bloc+1].a[bloc+1].coefsT
            hi = h[i+1]
            for j in 0:(N-1)
                target[j+1] += mu[j+1] * hi
            end
        end
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
function tGswEncryptZero(result::TGswSample, alpha::Float64, key::TGswKey)
    rlkey = key.tlwe_key
    kpl = key.params.kpl

    for p in 0:(kpl-1)
        tLweSymEncryptZero(result.all_sample[p+1], alpha, rlkey)
    end
end


# mult externe de X^{a_i} par bki
function tGswMulByXaiMinusOne(result::TGswSample, ai::Int32, bk::TGswSample, params::TGswParams)
    par = params.tlwe_params
    kpl = params.kpl
    for i in 0:(kpl-1)
        tLweMulByXaiMinusOne(result.all_sample[i+1], ai, bk.all_sample[i+1], par)
    end
end


# Update l'accumulateur ligne 5 de l'algo toujours
function tGswExternMulToTLwe(accum::TLweSample, sample::TGswSample, params::TGswParams)
    par = params.tlwe_params
    N = par.N
    kpl = params.kpl
    # TODO: improve this new/delete
    dec = new_IntPolynomial_array(kpl, N)

    tGswTLweDecompH(dec, accum, params)
    tLweClear(accum, par)
    for i in 0:(kpl-1)
        tLweAddMulRTo(accum, dec[i+1], sample.all_sample[i+1], par)
    end
end

# encrypts a poly message
function tGswSymEncrypt(result::TGswSample, message::IntPolynomial, alpha::Float64, key::TGswKey)
    tGswEncryptZero(result, alpha, key)
    tGswAddMuH(result, message, key.params)
end


# encrypts a constant message
function tGswSymEncryptInt(result::TGswSample, message::Int32, alpha::Float64, key::TGswKey)
    tGswEncryptZero(result, alpha, key)
    tGswAddMuIntH(result, message, key.params)
end


# encrypts a message = 0 ou 1
function tGswEncryptB(result::TGswSample, message::Int32, alpha::Float64, key::TGswKey)
    tGswEncryptZero(result, alpha, key)
    if (message == 1)
        tGswAddH(result, key.params)
    end
end


# à revoir
function tGswSymDecrypt(result::IntPolynomial, sample::TGswSample, key::TGswKey, Msize::Int32)
    params = key.params
    rlwe_params = params.tlwe_params
    N = rlwe_params.N
    l = params.l
    k = rlwe_params.k
    testvec = TorusPolynomial(N)
    tmp = TorusPolynomial(N)

    indic = modSwitchToTorus32(1, Msize)
    torusPolynomialClear(testvec)
    testvec.coefsT[0+1] = indic
    decomp = tGswTorus32PolynomialDecompH(testvec, params)

    torusPolynomialClear(testvec)
    for i in 0:(l-1)
        for j in 1:(N-1)
            @assert decomp[i+1].coefs[j+1] == 0
        end
        tLwePhase(tmp, sample.bloc_sample[i+1, k+1], key.tlwe_key)
        torusPolynomialAddMulR(testvec, decomp[i+1], tmp)
    end
    for i in 0:(N-1)
        result.coefs[i+1] = modSwitchFromTorus32(testvec.coefsT[i+1], Msize)
    end
end


# à revoir
function tGswSymDecryptInt(sample::TGswSample, key::TGswKey)
    phase = TorusPolynomial(key.params.tlwe_params.N)

    tGswPhase(phase, sample, key)
    modSwitchFromTorus32(phase.coefsT[0+1], Msize)
end


# fonction de decomposition
function tGswTLweDecompH(result::Array{IntPolynomial, 1}, sample::TLweSample, params::TGswParams)
    k = params.tlwe_params.k
    l = params.l

    for i in 0:k # b=a[k]
        result[i*l+1:(i+1)*l] = tGswTorus32PolynomialDecompH(sample.a[i], params)
    end
end


function Torus32PolynomialDecompH_old(result::IntPolynomial, sample::TorusPolynomial, params::TGswParams)
    N = params.tlwe_params.N
    l = params.l
    Bgbit = params.Bgbit
    maskMod = params.maskMod
    halfBg = params.halfBg
    offset = params.offset

    for j in 0:(N-1)
        temp0::UInt32 = sample.coefsT[j+1] + offset
        for p in 0:(l-1)
            temp1::UInt32 = (temp0 >> (32 - (p + 1) * Bgbit)) & maskMod; # doute
            result[p+1].coefs[j+1] = temp1 - halfBg
        end
    end
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


#result = a*b
function tGswExternProduct(result::TLweSample, a::TGswSample, b::TLweSample, params::TGswParams)
    parlwe = params.tlwe_params
    N = parlwe.N
    kpl = params.kpl
    dec = new_IntPolynomial_array(kpl, N)

    tGswTLweDecompH(dec, b, params)
    tLweClear(result, parlwe)

    for i in 0:(kpl-1)
        tLweAddMulRTo(result, dec[i+1], a.all_sample[i+1], parlwe)
    end

    result.current_variance += b.current_variance # todo + the error term?
end


# result = (0,mu)
function tGswNoiselessTrivial(result::TGswSample, mu::IntPolynomial, params::TGswParams)
    tGswClear(result, params)
    tGswAddMuH(result, mu, params)
end


# For all the kpl TLWE samples composing the TGSW sample
# It computes the inverse FFT of the coefficients of the TLWE sample
function tGswToFFTConvert(result::TGswSampleFFT, source::TGswSample, params::TGswParams)
    kpl = params.kpl

    for p in 0:(kpl-1)
        tLweToFFTConvert(result.all_samples[p+1], source.all_sample[p+1], params.tlwe_params)
    end
end

# For all the kpl TLWE samples composing the TGSW sample
# It computes the FFT of the coefficients of the TLWEfft sample
function tGswFromFFTConvert(result::TGswSample, source::TGswSampleFFT, params::TGswParams)
    kpl = params.kpl

    for p in 0:(kpl-1)
        tLweFromFFTConvert(result.all_sample[p+1], source.all_samples[p+1], params.tlwe_params)
    end
end


# result = result + H
function tGswFFTAddH(result::TGswSampleFFT, params::TGswParams)
    k = params.tlwe_params.k
    l = params.l

    for j in 0:(l-1)
        hj = params.h[j+1]
        for i in 0:k
            LagrangeHalfCPolynomialAddTorusConstant(result.sample[j+1,i+1].a[i+1], hj)
        end
    end
end

# result = list of TLWE (0,0)
function tGswFFTClear(result::TGswSampleFFT, params::TGswParams)
    kpl = params.kpl
    for p in 0:(kpl-1)
        tLweFFTClear(result.all_samples[p+1], params.tlwe_params)
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
