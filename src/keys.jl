struct TFHEParameters
    ks_decomp_length :: Int
    ks_log2_base :: Int
    in_out_params :: LweParams
    tgsw_params :: TGswParams

    function TFHEParameters()
        # the parameters are only implemented for about 128bit of security!

        lwe_size = 500

        tlwe_polynomial_degree = 1024
        tlwe_mask_size = 1

        bs_decomp_length = 2 # bootstrap decomposition length
        bs_log2_base = 10 # bootstrap log2(decomposition_base)

        ks_decomp_length = 8 # keyswitch decomposition length (the precision of the keyswitch)
        ks_log2_base = 2 # keyswitch log2(decomposition base)

        ks_stdev = 1/2^15 * sqrt(2 / pi) # keyswitch minimal standard deviation
        bs_stdev = 9e-9 * sqrt(2 / pi) # bootstrap minimal standard deviation
        max_stdev = 1/2^4 / 4 * sqrt(2 / pi) # max standard deviation for a 1/4 msg space

        params_in = LweParams(lwe_size, ks_stdev, max_stdev)
        params_accum = TLweParams(tlwe_polynomial_degree, tlwe_mask_size, bs_stdev, max_stdev)
        params_bs = TGswParams(bs_decomp_length, bs_log2_base, params_accum)

        new(ks_decomp_length, ks_log2_base, params_in, params_bs)
    end
end


struct TFHESecretKey
    params :: TFHEParameters
    lwe_key :: LweKey

    function TFHESecretKey(rng::AbstractRNG, params::TFHEParameters)
        lwe_key = LweKey(rng, params.in_out_params)
        new(params, lwe_key)
    end
end


struct TFHECloudKey
    params :: TFHEParameters
    bootstrap_key :: BootstrapKey
    keyswitch_key :: KeyswitchKey

    function TFHECloudKey(rng::AbstractRNG, secret_key::TFHESecretKey)
        params = secret_key.params
        tgsw_key = TGswKey(rng, params.tgsw_params)

        bs_key = BootstrapKey(rng, secret_key.lwe_key, tgsw_key)
        ks_key = KeyswitchKey(
            rng, params.ks_decomp_length, params.ks_log2_base, secret_key.lwe_key, tgsw_key)

        new(secret_key.params, bs_key, ks_key)
    end
end


tfhe_parameters(key::TFHESecretKey) = key.params
tfhe_parameters(key::TFHECloudKey) = key.params


function tfhe_key_pair(rng::AbstractRNG, params::Union{Nothing, TFHEParameters}=nothing)
    if params === nothing
        params = TFHEParameters()
    end
    secret_key = TFHESecretKey(rng, params)
    cloud_key = TFHECloudKey(rng, secret_key)
    secret_key, cloud_key
end


# encrypts a boolean
function tfhe_encrypt_bit(rng::AbstractRNG, key::TFHESecretKey, message::Bool)
    _1s8::Torus32 = modSwitchToTorus32(1, 8)
    mu::Torus32 = message ? _1s8 : -_1s8
    alpha = key.params.in_out_params.min_noise # TODO: specify noise
    lwe_encrypt(rng, mu, alpha, key.lwe_key)
end


# decrypts a boolean
function tfhe_decrypt_bit(key::TFHESecretKey, sample::TFHEEncryptedBit)
    lwe_phase(sample, key.lwe_key) > 0
end

TFHEEncryptedBit(params::TFHEParameters) = LweSample(params.in_out_params)
