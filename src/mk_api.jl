mktfhe_parameters_2party = SchemeParameters(
    500, 0.012467, # LWE parameters
    1024, 1, # TLWE parameters
    4, 7, 3.29e-10, # bootstrap parameters
    8, 2, 2.44e-5, # keyswitch parameters
    2
    )

mktfhe_parameters_4party = SchemeParameters(
    500, 0.012467, # LWE parameters
    1024, 1, # TLWE parameters
    5, 6, 3.29e-10, # bootstrap parameters
    8, 2, 2.44e-5, # keyswitch parameters
    4
    )

mktfhe_parameters_8party = SchemeParameters(
    500, 0.012467, # LWE parameters
    1024, 1, # TLWE parameters
    8, 4, 3.29e-10, # bootstrap parameters
    8, 2, 2.44e-5, # keyswitch parameters
    8
    )


# Resolving a circular dependency. SharedKey is used internally,
# and we don't want to expose SchemeParameters there.
function SharedKey(rng::AbstractRNG, params::SchemeParameters)
    tgsw_params = tgsw_parameters(params)
    tlwe_params = tlwe_parameters(params)
    SharedKey(rng, tgsw_params, tlwe_params)
end


# A part of the cloud key (bootstrap key + keyswitch key) generated independently by each party
# (since it involves their secret keys).
struct CloudKeyPart
    params :: SchemeParameters
    bk_part :: BootstrapKeyPart
    ks :: KeyswitchKey

    function CloudKeyPart(rng, secret_key::SecretKey, shared_key::SharedKey)
        params = secret_key.params
        tgsw_params = tgsw_parameters(params)
        tlwe_key = TLweKey(rng, tlwe_parameters(params))
        pk = PublicKey(rng, tlwe_key, params.bs_noise_stddev, shared_key, tgsw_params)
        bk = BootstrapKeyPart(rng, secret_key.key, tlwe_key, params.bs_noise_stddev, shared_key, pk)
        ks = KeyswitchKey(
                rng, params.ks_noise_stddev, keyswitch_parameters(params),
                secret_key.key, tlwe_key)
        new(params, bk, ks)
    end
end


struct MKCloudKey
    parties :: Int
    params :: SchemeParameters
    bootstrap_key :: MKBootstrapKey
    keyswitch_key :: Array{KeyswitchKey, 1}

    function MKCloudKey(ck_parts::Array{CloudKeyPart, 1})
        params = ck_parts[1].params
        parties = length(ck_parts)
        @assert parties <= params.max_parties

        bk = MKBootstrapKey([part.bk_part for part in ck_parts])
        ks = [part.ks for part in ck_parts]

        new(parties, params, bk, ks)
    end
end


function mk_encrypt(rng, secret_keys::Array{SecretKey, 1}, message::Bool)

    # TODO: (issue #6) encrypt separately for each party and share <a_i, s_i>?

    mu = encode_message(message ? 1 : -1, 8)

    params = secret_keys[1].params
    lwe_params = lwe_parameters(params)
    alpha = params.lwe_noise_stddev
    parties = length(secret_keys)

    a = hcat([rand_uniform_torus32(rng, lwe_params.size) for i in 1:parties]...)
    b = (rand_gaussian_torus32(rng, mu, alpha)
        + reduce(+, a .* hcat([secret_keys[i].key.key for i in 1:parties]...)))

    MKLweSample(lwe_params, a, b, alpha^2)
end


function mk_decrypt(secret_keys::Array{SecretKey, 1}, sample::MKLweSample)
    # TODO: (issue #6) decrypt separately at each party and join phases?
    mk_lwe_phase(sample, [sk.key for sk in secret_keys]) > 0
end
