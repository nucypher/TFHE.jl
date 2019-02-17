struct KeyswitchKey

    input_size :: Int32 # length of the input key: s'
    decomp_length :: Int32 # decomposition length
    log2_base :: Int32 # log_2(base)
    out_params :: LweParams # params of the output key s
    key :: Array{LweSample, 3} # the keyswitch elements: a n.l.base matrix

    function KeyswitchKey(
            rng::AbstractRNG, ks_decomp_length::Int, ks_log2_base::Int, key_in::LweKey, tgsw_key::TGswKey)

        bk_params = tgsw_key.params
        in_out_params = key_in.params
        accum_params = bk_params.tlwe_params
        extract_params = accum_params.extracted_lweparams

        n = in_out_params.len
        N = extract_params.len

        accum_key = tgsw_key.tlwe_key
        extracted_key = LweKey(extract_params, accum_key)

        in_key = extracted_key
        out_key = key_in
        n = N

        out_params = out_key.params
        t = ks_decomp_length
        basebit = ks_log2_base

        base::Int32 = 1 << basebit
        ks0_raw = [LweSample(out_params) for i in 1:(n * t * base)]
        ks = reshape(ks0_raw, Int64(base), Int64(t), Int64(n))

        alpha::Float64 = out_key.params.min_noise
        sizeks::Int32 = n * t * (base - 1)

        err::Float64 = 0

        # chose a random vector of gaussian noises
        noise = Array{Float64}(undef, sizeks)
        for i in 0:(sizeks-1)
            noise[i+1] = rand_gaussian_float(rng, alpha)
            err += noise[i+1]
        end
        # recenter the noises
        err = err / sizeks
        for i in 0:(sizeks-1)
            noise[i+1] -= err
        end

        # generate the ks
        index :: Int32 = 0
        for i in 0:(n-1)
            for j in 0:(t-1)

                # term h=0 as trivial encryption of 0 (it will not be used in the KeySwitching)
                ks[0+1,j+1,i+1] = lwe_noiseless_trivial(Torus32(0), out_key.params)
                for h in 1:(base-1) # pas le terme en 0
                    mess::Torus32 = (in_key.key[i+1] * Int32(h)) * Int32(1 << (32 - (j + 1) * basebit))
                    ks[h+1,j+1,i+1] = lwe_encrypt(
                        rng, mess, noise[index+1], alpha, out_key)
                    index += 1
                end
            end
        end

        new(n, t, basebit, out_params, ks)
    end
end



#=
 * translates the message of the result sample by -sum(a[i].s[i]) where s is the secret
 * embedded in ks.
 * @param result the LWE sample to translate by -sum(ai.si).
 * @param ks The (n x t x base) key switching key
 *        ks[i][j][k] encodes k.s[i]/base^(j+1)
 * @param params The common LWE parameters of ks and result
 * @param ai The input torus array
 * @param n The size of the input key
 * @param t The precision of the keyswitch (technically, 1/2.base^t)
 * @param basebit Log_2 of base
=#
function lweKeySwitchTranslate_fromArray(result::LweSample,
        ks::Array{LweSample, 3}, params::LweParams,
        ai::Array{Torus32, 1},
        n::Int32, t::Int32, basebit::Int32)

    base::Int32 = 1 << basebit # base=2 in [CGGI16]
    prec_offset::Int32 = 1 << (32 - (1 + basebit * t)) # precision
    mask::Int32 = base - 1

    for i in 0:(n-1)
        aibar::UInt32 = unsigned(ai[i+1] + prec_offset)
        for j in 0:(t-1)
            aij::UInt32 = (aibar >> (32 - (j + 1) * basebit)) & mask
            if aij != 0
                result -= ks[aij+1,j+1,i+1]
            end
        end
    end

    result
end


#sample=(a',b')
function lweKeySwitch(ks::KeyswitchKey, sample::LweSample)
    params = ks.out_params
    n = ks.input_size
    basebit = ks.log2_base
    t = ks.decomp_length

    result = lwe_noiseless_trivial(sample.b, params)
    lweKeySwitchTranslate_fromArray(result, ks.key, params, sample.a, n, t, basebit)
end
