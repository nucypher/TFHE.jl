struct _LweKeySwitchKey
    len :: Int32 # length of the input key: s'
    t :: Int32 # decomposition length
    basebit :: Int32 # log_2(base)
    base :: Int32 # decomposition base: a power of 2
    out_params :: LweParams # params of the output key s
    ks :: Array{LweSample, 3} # the keyswitch elements: a n.l.base matrix

    #=
    Create the key switching key:
     * normalize the error in the beginning
     * chose a random vector of gaussian noises (same size as ks)
     * recenter the noises
     * generate the ks by creating noiseless encryprions and then add the noise
    =#
    function _LweKeySwitchKey(
            rng::AbstractRNG,
            n::Int, t::Int, basebit::Int,
            in_key::LweKey, out_key::LweKey)

        out_params = out_key.params

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
                ks[0+1,j+1,i+1] = lweNoiselessTrivial(Torus32(0), out_key.params)
                for h in 1:(base-1) # pas le terme en 0
                    mess::Torus32 = (in_key.key[i+1] * Int32(h)) * Int32(1 << (32 - (j + 1) * basebit))
                    lweSymEncryptWithExternalNoise(
                        rng, ks[h+1,j+1,i+1], mess, noise[index+1], alpha, out_key)
                    index += 1
                end
            end
        end

        new(n, t, basebit, base, out_params, ks)
    end
end


struct KeyswitchKey

    len :: Int32 # length of the input key: s'
    t :: Int32 # decomposition length
    basebit :: Int32 # log_2(base)
    base :: Int32 # decomposition base: a power of 2
    out_params :: LweParams # params of the output key s
    # these don't seem to be used anywhere
    #ks0_raw :: Array{LweSample, 1} # tableau qui contient tout les Lwe samples de taille nlbase
    #ks1_raw :: Array{LweSample, 2} # de taille nl  pointe vers un tableau ks0_raw dont les cases sont espaceés de base positions
    ks :: Array{LweSample, 3} # the keyswitch elements: a n.l.base matrix

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

        ks = _LweKeySwitchKey(rng, N, ks_decomp_length, ks_log2_base, extracted_key, key_in)

        new(ks.len, ks.t, ks.basebit, ks.base, ks.out_params, ks.ks)
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
    n = ks.len
    basebit = ks.basebit
    t = ks.t

    result = lweNoiselessTrivial(sample.b, params)
    lweKeySwitchTranslate_fromArray(result, ks.ks, params, sample.a, n, t, basebit)
end
