struct BootstrapKey
    in_out_params :: LweParams # paramètre de l'input et de l'output. key: s
    extract_params :: LweParams # params after extraction: key: s'
    bk_params :: TGswParams # params of the Gsw elems in bk. key: s"
    key :: Array{TransformedTGswSample, 1} # the bootstrapping key (s->s")

    function BootstrapKey(rng::AbstractRNG, lwe_key::LweKey, tgsw_key::TGswKey)

        in_out_params = lwe_key.params
        bk_params = tgsw_key.params
        accum_params = bk_params.tlwe_params
        extract_params = accum_params.extracted_lweparams

        kin = lwe_key.key
        alpha = accum_params.min_noise
        lwe_len = in_out_params.len

        bk = [tgsw_encrypt(rng, kin[i], alpha, tgsw_key) for i in 1:lwe_len]
        transformed_bk = forward_transform.(bk)

        new(in_out_params, extract_params, bk_params, transformed_bk)
    end
end


function mux_rotate(
        accum::TLweSample, bki::TransformedTGswSample, barai::Int32, bk_params::TGswParams)
    # accum += BK_i * [(X^bar{a}_i-1) * accum]
    temp = shift_polynomial(accum, barai) - accum
    accum + tgsw_extern_mul(temp, bki, bk_params)
end


#=
 * multiply the accumulator by X^sum(bara_i.s_i)
 * @param accum the TLWE sample to multiply
 * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
 * @param bara An array of n coefficients between 0 and 2N-1
=#
function blind_rotate(accum::TLweSample, bk::BootstrapKey, bara::Array{Int32, 1})
    for i in 1:length(bk.key)
        if bara[i] != 0
            accum = mux_rotate(accum, bk.key[i], bara[i], bk.bk_params)
        end
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
=#
function blind_rotate_and_extract(
        v::TorusPolynomial, bk::BootstrapKey, barb::Int32, bara::Array{Int32, 1})

    accum_params = bk.bk_params.tlwe_params

    # testvector = X^{2N-barb}*v == X^{-barb}*v
    testvectbis = shift_polynomial(v, -barb)

    accum = tlwe_noiseless_trivial(testvectbis, accum_params)
    accum = blind_rotate(accum, bk, bara)
    tlwe_extract_sample(accum, accum_params)
end


#=
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
=#
function bootstrap_wo_keyswitch(bk::BootstrapKey, mu::Torus32, x::LweSample)

    p_degree = bk.bk_params.tlwe_params.polynomial_degree

    # Modulus switching
    bara = modSwitchFromTorus32.(x.a, p_degree * 2)
    barb = modSwitchFromTorus32(x.b, p_degree * 2)

    # the initial testvec = [mu,mu,mu,...,mu]
    testvect = torus_polynomial(repeat([mu], p_degree))

    # Bootstrapping rotation and extraction
    blind_rotate_and_extract(testvect, bk, barb, bara)
end


#=
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
=#
function bootstrap(bk::BootstrapKey, ks::KeyswitchKey, mu::Torus32, x::LweSample)
    u = bootstrap_wo_keyswitch(bk, mu, x)
    keyswitch(ks, u)
end
