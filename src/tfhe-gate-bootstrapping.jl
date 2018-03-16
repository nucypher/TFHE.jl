#=
tfhe_gate_bootstrapping_structures.h
tfhe_gate_bootstrapping_functions.h
tfhe_gate_bootstrapping_structures.cpp
thfe_gate_bootstrapping.cpp
=#

struct TFheGateBootstrappingParameterSet
    ks_t :: Int32
    ks_basebit :: Int32
    in_out_params :: LweParams
    tgsw_params :: TGswParams
end


struct TFheGateBootstrappingCloudKeySet
    params :: TFheGateBootstrappingParameterSet
    bk :: LweBootstrappingKey
    bkFFT :: LweBootstrappingKeyFFT
end


struct TFheGateBootstrappingSecretKeySet
    params :: TFheGateBootstrappingParameterSet
    lwe_key :: LweKey
    tgsw_key :: TGswKey
    cloud :: TFheGateBootstrappingCloudKeySet

    TFheGateBootstrappingSecretKeySet(
        params::TFheGateBootstrappingParameterSet,
        bk::LweBootstrappingKey,
        bkFFT::LweBootstrappingKeyFFT,
        lwe_key::LweKey,
        tgsw_key::TGswKey) = new(
            params, lwe_key, tgsw_key,
            TFheGateBootstrappingCloudKeySet(params, bk, bkFFT))
end



mulBySqrtTwoOverPi(x) = x * sqrt(2 / pi)


function new_default_gate_bootstrapping_parameters(minimum_lambda::Int)

    if minimum_lambda > 128
        error("Sorry, for now, the parameters are only implemented for about 128bit of security!")
    end

    N::Int32 = 1024
    k::Int32 = 1
    n::Int32 = 500
    bk_l::Int32 = 2
    bk_Bgbit::Int32 = 10
    ks_basebit::Int32 = 2
    ks_length::Int32 = 8
    ks_stdev = mulBySqrtTwoOverPi(1/2^15) # standard deviation
    bk_stdev = mulBySqrtTwoOverPi(9e-9) # standard deviation
    max_stdev = mulBySqrtTwoOverPi(1/2^4 / 4) # max standard deviation for a 1/4 msg space

    params_in = LweParams(n, ks_stdev, max_stdev)
    params_accum = TLweParams(N, k, bk_stdev, max_stdev)
    params_bk = TGswParams(bk_l, bk_Bgbit, params_accum)

    TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk)
end


function new_random_gate_bootstrapping_secret_keyset(params::TFheGateBootstrappingParameterSet)
    lwe_key = LweKey(params.in_out_params)
    lweKeyGen(lwe_key)
    tgsw_key = TGswKey(params.tgsw_params)
    tGswKeyGen(tgsw_key)
    bk = LweBootstrappingKey(
        params.ks_t, params.ks_basebit, params.in_out_params, params.tgsw_params)
    tfhe_createLweBootstrappingKey(bk, lwe_key, tgsw_key)
    bkFFT = LweBootstrappingKeyFFT(bk)

    TFheGateBootstrappingSecretKeySet(params, bk, bkFFT, lwe_key, tgsw_key)
end


# generate a new unititialized ciphertext
function new_gate_bootstrapping_ciphertext(params::TFheGateBootstrappingParameterSet)
    LweSample(params.in_out_params)
end

# generate a new unititialized ciphertext array of length nbelems
function new_gate_bootstrapping_ciphertext_array(nbelems::Int, params::TFheGateBootstrappingParameterSet)
    [LweSample(params.in_out_params) for i in 1:nbelems]
end

# encrypts a boolean
function bootsSymEncrypt(result::LweSample, message::Bool, key::TFheGateBootstrappingSecretKeySet)
    _1s8::Torus32 = modSwitchToTorus32(1, 8)
    mu::Torus32 = message ? _1s8 : -_1s8
    alpha = key.params.in_out_params.alpha_min # TODO: specify noise
    lweSymEncrypt(result, mu, alpha, key.lwe_key)
end

# decrypts a boolean
function bootsSymDecrypt(sample::LweSample, key::TFheGateBootstrappingSecretKeySet)
    mu = lwePhase(sample, key.lwe_key)
    (mu > 0 ? Int32(1) : Int32(0)) # we have to do that because of the C binding
end
