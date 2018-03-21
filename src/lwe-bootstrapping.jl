function lwe_bootstrapping_key(
        rng::AbstractRNG, ks_t::Int32, ks_basebit::Int32, key_in::LweKey, rgsw_key::TGswKey)

    bk_params = rgsw_key.params
    in_out_params = key_in.params
    accum_params = bk_params.tlwe_params
    extract_params = accum_params.extracted_lweparams

    n = in_out_params.n
    N = extract_params.n

    accum_key = rgsw_key.tlwe_key
    extracted_key = LweKey(extract_params, accum_key)

    ks = LweKeySwitchKey(rng, N, ks_t, ks_basebit, extracted_key, key_in)

    bk = TGswSampleArray(bk_params, Int(n))
    kin = key_in.key
    alpha = accum_params.alpha_min

    rng2 = MersenneTwister(123) # !!!
    tGswSymEncryptInt(rng, rng2, bk, kin, alpha, rgsw_key)

    bk, ks
end


struct LweBootstrappingKeyFFT
    in_out_params :: LweParams # paramètre de l'input et de l'output. key: s
    bk_params :: TGswParams # params of the Gsw elems in bk. key: s"
    accum_params :: TLweParams # params of the accum variable key: s"
    extract_params :: LweParams # params after extraction: key: s'
    bkFFT :: TGswSampleFFTArray # the bootstrapping key (s->s")
    ks :: LweKeySwitchKey # the keyswitch key (s'->s)

    function LweBootstrappingKeyFFT(
            rng::AbstractRNG, ks_t::Int32, ks_basebit::Int32, lwe_key::LweKey, tgsw_key::TGswKey)

        in_out_params = lwe_key.params
        bk_params = tgsw_key.params
        accum_params = bk_params.tlwe_params
        extract_params = accum_params.extracted_lweparams

        bk, ks = lwe_bootstrapping_key(rng, ks_t, ks_basebit, lwe_key, tgsw_key)

        n = in_out_params.n

        # Bootstrapping Key FFT
        bkFFT = TGswSampleFFTArray(bk_params, Int(n))
        tGswToFFTConvert(bkFFT, bk, bk_params)

        new(in_out_params, bk_params, accum_params, extract_params, bkFFT, ks)
    end

end


function tfhe_MuxRotate_FFT(
        result::TLweSampleArray, accum::TLweSampleArray, bki::TGswSampleFFTArray, barai::Array{Int32},
        bk_params::TGswParams)

    # ACC = BKi*[(X^barai-1)*ACC]+ACC
    # temp = (X^barai-1)*ACC

    tLweMulByXaiMinusOne(result, barai, accum, bk_params.tlwe_params)

    # temp *= BKi
    tGswFFTExternMulToTLwe(result, bki, bk_params)

    # ACC += temp
    tLweAddTo(result, accum, bk_params.tlwe_params)
end


#=
 * multiply the accumulator by X^sum(bara_i.s_i)
 * @param accum the TLWE sample to multiply
 * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
=#
function tfhe_blindRotate_FFT(accum::TLweSampleArray,
                                 bkFFT::TGswSampleFFTArray,
                                 bara::Array{Int32},
                                 n::Int32,
                                 bk_params::TGswParams)

    #TGswSampleFFT* temp = new_TGswSampleFFT(bk_params);
    temp = TLweSampleArray(bk_params.tlwe_params, size(accum)...)
    temp2 = temp
    temp3 = accum

    accum_in_temp3 = true

    for i in 0:(n-1)

        barai = bara[i+1,:] # !!! assuming the ciphertext is 1D

        # TODO: is it possible to avoid work for those cases and still keep array arithmetic?
        #       Does it even save that much work?
        #if barai == 0
        #    continue #indeed, this is an easy case!
        #end

        tfhe_MuxRotate_FFT(temp2, temp3, view(bkFFT, i+1), barai, bk_params)

        temp2, temp3 = temp3, temp2
        accum_in_temp3 = !accum_in_temp3
    end

    if !accum_in_temp3 # temp3 != accum
        tLweCopy(accum, temp3, bk_params.tlwe_params)
    end
end


#=
 * result = LWE(v_p) where p=barb-sum(bara_i.s_i) mod 2N
 * @param result the output LWE sample
 * @param v a 2N-elt anticyclic function (represented by a TorusPolynomial)
 * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
 * @param barb A coefficients between 0 and 2N-1
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
=#
function tfhe_blindRotateAndExtract_FFT(result::LweSampleArray,
                                           v::TorusPolynomialArray,
                                           bk::TGswSampleFFTArray,
                                           barb::Array{Int32},
                                           bara::Array{Int32},
                                           n::Int32,
                                           bk_params::TGswParams)

    accum_params = bk_params.tlwe_params
    extract_params = accum_params.extracted_lweparams
    N = accum_params.N
    _2N = Int32(2) * N

    # Test polynomial
    testvectbis = TorusPolynomialArray(Int(N), size(result)...) # TODO: get rid of Int()
    # Accumulator
    acc = TLweSampleArray(accum_params, size(result)...)

    # testvector = X^{2N-barb}*v
    tp_mul_by_xai!(testvectbis, _2N - barb, v)

    tLweNoiselessTrivial(acc, testvectbis, accum_params)

    # Blind rotation
    tfhe_blindRotate_FFT(acc, bk, bara, n, bk_params)

    # Extraction
    tLweExtractLweSample(result, acc, extract_params, accum_params)
end


#=
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
=#
function tfhe_bootstrap_woKS_FFT(result::LweSampleArray,
                                    bk::LweBootstrappingKeyFFT,
                                    mu::Torus32,
                                    x::LweSampleArray)

    bk_params = bk.bk_params
    accum_params = bk.accum_params
    in_params = bk.in_out_params
    N = accum_params.N
    Nx2 = 2 * N
    n = in_params.n

    testvect = TorusPolynomialArray(Int(N), size(result)...) # TODO: get rid of Int()

    # Modulus switching
    barb = modSwitchFromTorus32.(x.b, Nx2)
    bara = modSwitchFromTorus32.(x.a, Nx2)

    # the initial testvec = [mu,mu,mu,...,mu]
    # TODO: use an appropriate method
    testvect.coefsT .= mu

    # Bootstrapping rotation and extraction
    tfhe_blindRotateAndExtract_FFT(result, testvect, bk.bkFFT, barb, bara, n, bk_params)
end


#=
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
=#
function tfhe_bootstrap_FFT(result::LweSampleArray,
                               bk::LweBootstrappingKeyFFT,
                               mu::Torus32,
                               x::LweSampleArray)

    u = LweSampleArray(bk.accum_params.extracted_lweparams, size(result)...)

    tfhe_bootstrap_woKS_FFT(u, bk, mu, x)

    # Key switching
    lweKeySwitch(result, bk.ks, u)
end
