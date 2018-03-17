struct TFHEParameters
    ks_t :: Int32
    ks_basebit :: Int32
    in_out_params :: LweParams
    tgsw_params :: TGswParams

    function TFHEParameters()
        # In the reference implementation there was a parameter `minimum_lambda` here,
        # which was unused.

        # the parameters are only implemented for about 128bit of security!

        mul_by_sqrt_two_over_pi(x) = x * sqrt(2 / pi)

        N::Int32 = 1024
        k::Int32 = 1
        n::Int32 = 500
        bk_l::Int32 = 2
        bk_Bgbit::Int32 = 10
        ks_basebit::Int32 = 2
        ks_length::Int32 = 8
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
    bk :: LweBootstrappingKey
    bkFFT :: LweBootstrappingKeyFFT
end


tfhe_parameters(key::TFHESecretKey) = key.params
tfhe_parameters(key::TFHECloudKey) = key.params


function tfhe_key_pair(rng::AbstractRNG)
    params = TFHEParameters()

    lwe_key = LweKey(params.in_out_params)
    lweKeyGen(lwe_key)
    tgsw_key = TGswKey(params.tgsw_params)
    tGswKeyGen(tgsw_key)
    secret_key = TFHESecretKey(params, lwe_key, tgsw_key)

    bk = LweBootstrappingKey(
        params.ks_t, params.ks_basebit, params.in_out_params, params.tgsw_params)
    tfhe_createLweBootstrappingKey(bk, lwe_key, tgsw_key)
    bkFFT = LweBootstrappingKeyFFT(bk)

    cloud_key = TFHECloudKey(params, bk, bkFFT)
    secret_key, cloud_key
end


# encrypts a boolean
function tfhe_encrypt_bit!(key::TFHESecretKey, result::TFHEEncryptedBit, message::Bool)
    _1s8::Torus32 = modSwitchToTorus32(1, 8)
    mu::Torus32 = message ? _1s8 : -_1s8
    alpha = key.params.in_out_params.alpha_min # TODO: specify noise
    lweSymEncrypt(result, mu, alpha, key.lwe_key)
end


# decrypts a boolean
function tfhe_decrypt_bit!(key::TFHESecretKey, sample::TFHEEncryptedBit)
    mu = lwePhase(sample, key.lwe_key)
    (mu > 0 ? Int32(1) : Int32(0))
end

TFHEEncryptedBit(params::TFHEParameters) = LweSample(params.in_out_params)
