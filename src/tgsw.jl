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


function tGswTorus32PolynomialDecompH(
        result::IntPolynomialArray, sample::TorusPolynomialArray, params::TGswParams)

    N = params.tlwe_params.N
    l = params.l
    k = params.tlwe_params.k
    Bgbit = params.Bgbit

    maskMod = params.maskMod
    halfBg = params.halfBg
    offset = params.offset

    decal(p) = 32 - p * Bgbit

    ps = reshape(1:l, 1, l, 1, 1)
    sample_coefs = reshape(sample.coefsT, N, 1, size(sample)...)

    # do the decomposition
    @. result.coefs = ((sample_coefs + offset) >> decal(ps)) & maskMod - halfBg
end


# For all the kpl TLWE samples composing the TGSW sample
# It computes the inverse FFT of the coefficients of the TLWE sample
function tGswToFFTConvert(result::TGswSampleFFTArray, source::TGswSampleArray, params::TGswParams)
    tLweToFFTConvert(result.samples, source.samples, params.tlwe_params)
end


function tLweFFTAddMulRTo(res, a, b, bk_idx)
    #=
    # FIXME: When Julia is smart enough, can be replaced by
    d = reshape(a, div(N, 2), 1, l, k+1, size(accum, 1))
    @views for i in 1:(k+1)
        for j in 1:l
            tmpa.a.coefsC .+= d[:,:,j,i,:] .* b[:,:,j,i,bk_idx]
        end
    end
    =#

    N2 = size(res, 1)
    k = size(res, 2)
    ml = size(res, 3)
    l = size(a, 2)
    @inbounds @simd for m in 1:ml
    for i in 1:k
        for j in 1:l
                for q in 1:k
                for p in 1:N2
                        res[p,q,m] = res[p,q,m] + a[p,j,i,m] * b[p,q,j,i,bk_idx]
                    end
                end
            end
        end
    end
end


# External product (*): accum = gsw (*) accum
function tGswFFTExternMulToTLwe(
        accum::TLweSampleArray, gsw::TGswSampleFFTArray, bk_idx, params::TGswParams,
        tmpa::TLweSampleFFTArray, deca::IntPolynomialArray, decaFFT::LagrangeHalfCPolynomialArray)

    tlwe_params = params.tlwe_params
    k = tlwe_params.k
    l = params.l
    kpl = params.kpl
    N = tlwe_params.N

    tGswTorus32PolynomialDecompH(deca, accum.a, params)

    ip_ifft!(decaFFT, deca)

    tLweFFTClear(tmpa, tlwe_params)

    N2 = div(N, 2)
    ml = size(accum, 1)
    res = tmpa.a.coefsC
    a = decaFFT.coefsC
    b = gsw.samples.a.coefsC

    tLweFFTAddMulRTo(res, a, b, bk_idx)

    tLweFromFFTConvert(accum, tmpa, tlwe_params)
end
