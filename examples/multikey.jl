using Random
using TFHE


function main()
    parties = 2

    params = mktfhe_parameters_2party

    rng = MersenneTwister()

    # Processed on clients' machines
    secret_keys = [SecretKey(rng, params) for i in 1:parties]

    # Created by the server
    shared_key = SharedKey(rng, params)

    # Processed on clients' machines
    ck_parts = [CloudKeyPart(rng, secret_key, shared_key) for secret_key in secret_keys]

    # Processed on the server.
    # `ck_parts` only contain information `public_keys`, `secret_keys` remain secret.
    cloud_key = MKCloudKey(ck_parts)

    for trial = 1:10

        mess1 = rand(Bool)
        mess2 = rand(Bool)
        out = !(mess1 && mess2)

        enc_mess1 = mk_encrypt(rng, secret_keys, mess1)
        enc_mess2 = mk_encrypt(rng, secret_keys, mess2)

        dec_mess1 = mk_decrypt(secret_keys, enc_mess1)
        dec_mess2 = mk_decrypt(secret_keys, enc_mess2)
        @assert mess1 == dec_mess1
        @assert mess2 == dec_mess2

        enc_out = mk_gate_nand(cloud_key, enc_mess1, enc_mess2)

        dec_out = mk_decrypt(secret_keys, enc_out)
        @assert out == dec_out

        println("Trial $trial: $mess1 NAND $mess2 = $out")
    end
end


main()
