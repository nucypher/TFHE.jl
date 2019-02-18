struct BootstrapKey
    in_out_params :: LweParams # paramètre de l'input et de l'output. key: s
    bk_params :: TGswParams # params of the Gsw elems in bk. key: s"
    accum_params :: TLweParams # params of the accum variable key: s"
    extract_params :: LweParams # params after extraction: key: s'
    key :: Array{TransformedTGswSample, 1} # the bootstrapping key (s->s")

    function BootstrapKey(rng::AbstractRNG, lwe_key::LweKey, tgsw_key::TGswKey)

        in_out_params = lwe_key.params
        bk_params = tgsw_key.params
        accum_params = bk_params.tlwe_params
        extract_params = accum_params.extracted_lweparams

        kin = lwe_key.key
        alpha = accum_params.min_noise
        n = in_out_params.len
        bk = [tgsw_encrypt(rng, kin[i], alpha, tgsw_key) for i in 1:n]

        # Bootstrapping Key FFT
        bkFFT = forward_transform.(bk)

        new(in_out_params, bk_params, accum_params, extract_params, bkFFT)
    end
end


function tfhe_MuxRotate_FFT(
        accum::TLweSample, bki::TransformedTGswSample, barai::Int32,
        bk_params::TGswParams)

    # ACC = BKi*[(X^barai-1)*ACC]+ACC

    # temp = (X^barai-1)*ACC
    result = shift_polynomial(accum, barai) - accum

    # temp *= BKi
    result = tgsw_extern_mul(result, bki, bk_params)

    # ACC += temp
    result += accum

    result
end


#=
 * multiply the accumulator by X^sum(bara_i.s_i)
 * @param accum the TLWE sample to multiply
 * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
=#
function tfhe_blindRotate_FFT(accum::TLweSample,
                                 bkFFT::Array{TransformedTGswSample, 1},
                                 bara::Array{Int32, 1},
                                 n::Int,
                                 bk_params::TGswParams)

    for i in 0:(n-1)

        barai = bara[i+1]
        if barai == 0
            continue
        end

        accum = tfhe_MuxRotate_FFT(accum, bkFFT[i+1], barai, bk_params)
    end

    accum
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
function tfhe_blindRotateAndExtract_FFT(
                                           v::TorusPolynomial,
                                           bk::Array{TransformedTGswSample, 1},
                                           barb::Int32,
                                           bara::Array{Int32, 1},
                                           n::Int,
                                           bk_params::TGswParams)

    accum_params = bk_params.tlwe_params

    # testvector = X^{2N-barb}*v == X^{-barb}*v
    if barb != 0
        testvectbis = shift_polynomial(v, -barb)
    else
        testvectbis = deepcopy(v)
    end

    # Accumulator
    acc = tlwe_noiseless_trivial(testvectbis, accum_params)

    # Blind rotation
    acc = tfhe_blindRotate_FFT(acc, bk, bara, n, bk_params)

    # Extraction
    tlwe_extract_sample(acc, accum_params)
end


#=
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
=#
function tfhe_bootstrap_woKS_FFT(bk::BootstrapKey, mu::Torus32, x::LweSample)

    bk_params = bk.bk_params
    accum_params = bk.accum_params
    in_params = bk.in_out_params
    p_degree = accum_params.polynomial_degree
    lwe_len = in_params.len

    bara = Array{Int32}(undef, lwe_len)

    # Modulus switching
    barb = modSwitchFromTorus32(x.b, p_degree * 2)
    for i in 0:(lwe_len-1)
        bara[i+1] = modSwitchFromTorus32(x.a[i+1], p_degree * 2)
    end

    # the initial testvec = [mu,mu,mu,...,mu]
    testvect = torus_polynomial(repeat([mu], p_degree))

    # Bootstrapping rotation and extraction
    tfhe_blindRotateAndExtract_FFT(testvect, bk.key, barb, bara, lwe_len, bk_params)
end


#=
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
=#
function tfhe_bootstrap_FFT(
                               bk::BootstrapKey,
                               ks::KeyswitchKey,
                               mu::Torus32,
                               x::LweSample)

    u = tfhe_bootstrap_woKS_FFT(bk, mu, x)
    lweKeySwitch(ks, u)
end
