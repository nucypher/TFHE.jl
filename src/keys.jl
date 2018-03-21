struct TFHEParameters
    ks_t :: Int
    ks_basebit :: Int
    in_out_params :: LweParams
    tgsw_params :: TGswParams

    function TFHEParameters()
        # In the reference implementation there was a parameter `minimum_lambda` here,
        # which was unused.

        # the parameters are only implemented for about 128bit of security!

        mul_by_sqrt_two_over_pi(x) = x * sqrt(2 / pi)

        N = 1024
        k = 1
        n = 500
        bk_l = 2
        bk_Bgbit = 10
        ks_basebit = 2
        ks_length = 8
        ks_stdev = mul_by_sqrt_two_over_pi(1/2^15) # standard deviation
        bk_stdev = mul_by_sqrt_two_over_pi(9e-9) # standard deviation
        max_stdev = mul_by_sqrt_two_over_pi(1/2^4 / 4) # max standard deviation for a 1/4 msg space

        params_in = LweParams(n, ks_stdev, max_stdev)
        params_accum = TLweParams(N, k, bk_stdev, max_stdev)
        params_bk = TGswParams(bk_l, bk_Bgbit, params_accum)

        new(ks_length, ks_basebit, params_in, params_bk)
    end
end


struct TFHESecretKey
    params :: TFHEParameters
    lwe_key :: LweKey
    tgsw_key :: TGswKey
end


struct TFHECloudKey
    params :: TFHEParameters
    bkFFT :: LweBootstrappingKeyFFT
end


tfhe_parameters(key::TFHESecretKey) = key.params
tfhe_parameters(key::TFHECloudKey) = key.params


function tfhe_key_pair(rng::AbstractRNG)
    params = TFHEParameters()

    lwe_key = LweKey(rng, params.in_out_params)
    tgsw_key = TGswKey(rng, params.tgsw_params)
    secret_key = TFHESecretKey(params, lwe_key, tgsw_key)

    bkFFT = LweBootstrappingKeyFFT(rng, params.ks_t, params.ks_basebit, lwe_key, tgsw_key)
    cloud_key = TFHECloudKey(params, bkFFT)

    secret_key, cloud_key
end


function tfhe_encrypt(rng::AbstractRNG, key::TFHESecretKey, message::BitArray)
    result = empty_ciphertext(key.params, length(message))
    _1s8 = modSwitchToTorus32(1, 8)
    mus = [bit ? _1s8 : -_1s8 for bit in message]
    alpha = key.params.in_out_params.alpha_min # TODO: specify noise
    lweSymEncrypt(rng, result, mus, alpha, key.lwe_key)
    result
end


function tfhe_decrypt(key::TFHESecretKey, ciphertext::LweSampleArray)
    mus = lwePhase(ciphertext, key.lwe_key)
    BitArray([(mu > 0) for mu in mus])
end


function empty_ciphertext(params::TFHEParameters, dims...)
    LweSampleArray(params.in_out_params, dims...)
end

