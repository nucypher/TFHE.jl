# Based on H. Chen, I. Chillotti, and Y. Song, Y"Multi-Key Homomophic Encryption from TFHE"

using Random
using TFHE:
    SchemeParameters, lwe_parameters, tlwe_parameters, tgsw_parameters, keyswitch_parameters,
    SecretKey, TLweKey,
    Torus32, IntPolynomial, TorusPolynomial, int_polynomial, torus_polynomial, transformed_mul,
    TransformedTorusPolynomial, inverse_transform,
    encode_message, decode_message,
    rand_uniform_bool, rand_uniform_torus32, rand_gaussian_float, dtot32, rand_gaussian_torus32,
    LweParams, TLweParams,
    LweKey, TLweKey, extract_lwe_key, lwe_phase,
    LweSample, KeyswitchKey, keyswitch, TGswParams, reverse_polynomial, decompose
import TFHE: forward_transform

using DarkIntegers


new_torus_polynomial(len) = TorusPolynomial(Array{Torus32}(undef, len), true)


mktfhe_parameters_2party = SchemeParameters(
    500, 0.012467, # LWE parameters
    1024, 1, # TLWE parameters
    4, 7, 3.29e-10, # bootstrap parameters
    8, 2, 2.44e-5, # keyswitch parameters
    2
    )


struct SharedKey

    params :: SchemeParameters
    a :: Array{TorusPolynomial, 1}

    function SharedKey(rng, params::SchemeParameters)
        decomp_length = tgsw_parameters(params).decomp_length
        p_degree = tlwe_parameters(params).polynomial_degree
        new(
            params,
            [torus_polynomial(rand_uniform_torus32(rng, p_degree)) for i in 1:decomp_length])
    end
end


struct PublicKey

    params :: SchemeParameters
    b :: Array{TorusPolynomial, 1}

    function PublicKey(rng, sk::SecretKey, tlwe_key::TLweKey, shared::SharedKey)

        params = sk.params
        decomp_length = tgsw_parameters(params).decomp_length
        p_degree = tlwe_parameters(params).polynomial_degree

        # The right part of the public key is `b_i = e_i + a*s_i`,
        # where `a` is shared between the parties
        # TODO: [1] was omitted in the original! It works while k=1, but will fail otherwise
        # TODO: this is basically tgsw_encrypt_zero() for mask_size=1
        b = [(
                transformed_mul(tlwe_key.key[1], shared.a[i])
                + torus_polynomial(
                    rand_gaussian_torus32(rng, zero(Int32), params.bs_noise_stddev, p_degree)))
            for i in 1:decomp_length]

        new(sk.params, b)
    end
end


# TGSW sample expanded for all the parties
# Since the C matrix returned by RGSW.Expand in the paper is sparse,
# we are keeping only nonzero elements of it.
mutable struct MKTGswExpSample

    tgsw_params :: TGswParams
    tlwe_params :: TLweParams

    x :: Array{TorusPolynomial, 2}
    y :: Array{TorusPolynomial, 2}
    c0 :: Array{TorusPolynomial, 1}
    c1 :: Array{TorusPolynomial, 1}

    current_variance :: Float64 # avg variance of the sample

    function MKTGswExpSample(tgsw_params::TGswParams, tlwe_params::TLweParams, parties::Int)

        decomp_length = tgsw_params.decomp_length
        p_degree = tlwe_params.polynomial_degree

        x = [new_torus_polynomial(p_degree) for i in 1:decomp_length, j in 1:parties]
        y = [new_torus_polynomial(p_degree) for i in 1:decomp_length, j in 1:parties]
        c0 = [new_torus_polynomial(p_degree) for i in 1:decomp_length]
        c1 = [new_torus_polynomial(p_degree) for i in 1:decomp_length]

        current_variance = 0.0

        new(tgsw_params, tlwe_params, x, y, c0, c1, current_variance)
    end
end


mutable struct MKTransformedTGswExpSample

    tgsw_params :: TGswParams
    tlwe_params :: TLweParams

    x :: Array{TransformedTorusPolynomial, 2}
    y :: Array{TransformedTorusPolynomial, 2}
    c0 :: Array{TransformedTorusPolynomial, 1}
    c1 :: Array{TransformedTorusPolynomial, 1}

    current_variance :: Float64 # avg variance of the sample

    MKTransformedTGswExpSample(tgsw_params, tlwe_params, x, y, c0, c1, current_variance) =
        new(tgsw_params, tlwe_params, x, y, c0, c1, current_variance)
end


function forward_transform(sample::MKTGswExpSample)
    MKTransformedTGswExpSample(
        sample.tgsw_params,
        sample.tlwe_params,
        forward_transform.(sample.x),
        forward_transform.(sample.y),
        forward_transform.(sample.c0),
        forward_transform.(sample.c1),
        sample.current_variance)
end


# Uni-encrypted TGSW sample
mutable struct MKTGswUESample

    tgsw_params :: TGswParams
    tlwe_params :: TLweParams

    c0 :: Array{TorusPolynomial, 1}
    c1 :: Array{TorusPolynomial, 1}
    d0 :: Array{TorusPolynomial, 1}
    d1 :: Array{TorusPolynomial, 1}
    f0 :: Array{TorusPolynomial, 1}
    f1 :: Array{TorusPolynomial, 1}

    current_variance :: Float64 # avg variance of the sample

    function MKTGswUESample(tgsw_params::TGswParams, tlwe_params::TLweParams)
        decomp_length = tgsw_params.decomp_length
        p_degree = tlwe_params.polynomial_degree

        c0 = [new_torus_polynomial(p_degree) for i in 1:decomp_length]
        c1 = [new_torus_polynomial(p_degree) for i in 1:decomp_length]
        d0 = [new_torus_polynomial(p_degree) for i in 1:decomp_length]
        d1 = [new_torus_polynomial(p_degree) for i in 1:decomp_length]
        f0 = [new_torus_polynomial(p_degree) for i in 1:decomp_length]
        f1 = [new_torus_polynomial(p_degree) for i in 1:decomp_length]

        current_variance = 0.0

        new(tgsw_params, tlwe_params, c0, c1, d0, d1, f0, f1, current_variance)
    end
end


# Encrypt an integer value
# In the paper: RGSW.UniEnc
# Similar to tgsw_encrypt()/tlwe_encrypt(), except the public key is supplied externally.
function mk_tgsw_encrypt(
        rng, message::Int32, alpha::Float64,
        secret_key::SecretKey, tlwe_key::TLweKey, shared_key::SharedKey, public_key::PublicKey)

    params = secret_key.params
    tgsw_params = tgsw_parameters(params)
    tlwe_params = tlwe_parameters(params)
    p_degree = tlwe_params.polynomial_degree
    decomp_length = tgsw_params.decomp_length

    result = MKTGswUESample(tgsw_params, tlwe_params)

    # The shared randomness
    r = int_polynomial(rand_uniform_bool(rng, p_degree))

    # C = (c0,c1) \in T^2dg, with c0 = s_party*c1 + e_c + m*g
    for i in 1:decomp_length
        result.c1[i] = torus_polynomial(rand_uniform_torus32(rng, p_degree))

        # c0 = s_party*c1 + e_c + m*g
        # TODO: it was just key, not key[1] in the original. Hardcoded mask_size=1? Check!
        result.c0[i] = (
            torus_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, p_degree))
            + message * tgsw_params.gadget_values[i]
            + transformed_mul(tlwe_key.key[1], result.c1[i]))
    end


    # D = (d0, d1) = r*[Pkey_party | Pkey_parties] + [E0|E1] + [0|m*g] \in T^2dg
    for i in 1:decomp_length

        # d1 = r*Pkey_parties[i] + E1 + m*g[i]
        result.d1[i] = (
            torus_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, p_degree))
            + message * tgsw_params.gadget_values[i]
            + transformed_mul(r, shared_key.a[i]))

        # d0 = r*Pkey_party[i] + E0
        result.d0[i] = (
            torus_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, p_degree))
            + transformed_mul(r, public_key.b[i]))
    end


    # F = (f0,f1) \in T^2dg, with f0 = s_party*f1 + e_f + r*g
    for i in 1:decomp_length

        result.f1[i] = torus_polynomial(rand_uniform_torus32(rng, p_degree))

        # f0 = s_party*f1 + e_f + r*g
        # TODO: it was just key, not key[1] in the original. Hardcoded mask_size=1? Check!
        result.f0[i] = (
            torus_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, p_degree))
            + r * tgsw_params.gadget_values[i]
            + transformed_mul(tlwe_key.key[1], result.f1[i]))
    end

    result.current_variance = alpha^2

    result
end


# In the paper: RGSW.Expand
function mk_tgsw_expand(sample::MKTGswUESample, party::Int, public_keys::Array{PublicKey, 1})

    tlwe_params = sample.tlwe_params
    tgsw_params = sample.tgsw_params
    p_degree = sample.tlwe_params.polynomial_degree
    decomp_length = tgsw_params.decomp_length
    parties = length(public_keys)

    result = MKTGswExpSample(tgsw_params, tlwe_params, parties)

    # Initialize: C = (0, ..., d1, ..., 0, c1, d0, ..., d0, ..., d0, c0)
    for j in 1:decomp_length
        for i in 1:parties
            result.y[j,i].coeffs .= 0
        end

        result.y[j,party].coeffs .= sample.d1[j].coeffs

        result.c1[j].coeffs .= sample.c1[j].coeffs

        for i in 1:parties
            result.x[j,i].coeffs .= sample.d0[j].coeffs
        end

        result.c0[j].coeffs .= sample.c0[j].coeffs
    end

    X = new_torus_polynomial(p_degree)
    Y = new_torus_polynomial(p_degree)
    b_temp = new_torus_polynomial(p_degree)
    u = [new_torus_polynomial(p_degree) for i in 1:decomp_length]

    for i in 1:parties
        if i != party
            for j in 1:decomp_length
                # b_temp = b_i[j] - b_party[j]
                b_temp = public_keys[i].b[j] - public_keys[party].b[j]
                # g^{-1}(b_temp) = [u_0, ...,u_dg-1]
                u = decompose(b_temp, tgsw_params)

                X.coeffs .= 0
                Y.coeffs .= 0
                for l in 1:decomp_length
                    # X = xi[j] = <g^{-1}(b_temp), f0>
                    X += transformed_mul(u[l], sample.f0[l])
                    # Y = yi[j] = <g^{-1}(b_temp), f1>
                    Y += transformed_mul(u[l], sample.f1[l])
                end

                # xi = d0 + xi
                result.x[j,i] += X
                # yi = 0 + yi
                result.y[j,i] += Y
            end
        end
    end

    # TODO: calculate the current variance correctly
    result.current_variance = sample.current_variance

    result
end


mutable struct MKLweSample

    params :: LweParams
    a :: Array{Torus32, 2} # masks from all parties: (n, parties)
    b :: Torus32 # the joined phase
    current_variance :: Float64 # average noise of the sample

    MKLweSample(params::LweParams, a::Array{Torus32, 2}, b::Torus32, cv::Float64) =
        new(params, a, b, cv)

    MKLweSample(params::LweParams, parties) =
        new(params, Array{Torus32}(undef, params.size, parties), Torus32(0), 0.)
end


Base.:-(x::MKLweSample, y::MKLweSample) =
    MKLweSample(x.params, x.a .- y.a, x.b - y.b, x.current_variance + y.current_variance)


Base.:+(x::MKLweSample, y::MKLweSample) =
    MKLweSample(x.params, x.a .+ y.a, x.b + y.b, x.current_variance + y.current_variance)


function mk_keyswitch(ks::Array{KeyswitchKey, 1}, sample::MKLweSample)
    parties = length(ks)
    result_parts = [
        keyswitch(ks[p], LweSample(sample.params, sample.a[:,p], Int32(0), 0.))
        for p in 1:parties]

    out_lwe_params = ks[1].out_lwe_params

    result = mk_lwe_noiseless_trivial(sample.b, out_lwe_params, parties)
    result + MKLweSample(
        out_lwe_params,
        hcat([part.a for part in result_parts]...),
        reduce(+, [part.b for part in result_parts]),
        0.0)
end


# A part of the bootstrap key generated independently by each party
# (since it involves their secret keys).
struct BootstrapKeyPart

    tgsw_params :: TGswParams
    tlwe_params :: TLweParams
    key_uni_enc :: Array{MKTGswUESample, 1}
    public_key :: PublicKey

    function BootstrapKeyPart(
            rng, secret_key::SecretKey, tlwe_key::TLweKey,
            shared_key::SharedKey, public_key::PublicKey)

        params = secret_key.params
        tgsw_params = tgsw_parameters(params)
        tlwe_params = tlwe_parameters(params)
        in_out_params = lwe_parameters(params)

        new(
            tgsw_params,
            tlwe_params,
            [mk_tgsw_encrypt(
                rng, secret_key.key.key[j],
                params.bs_noise_stddev, secret_key, tlwe_key, shared_key, public_key)
                for j in 1:in_out_params.size],
            public_key)
    end
end


# A part of the cloud key (bootstrap key + keyswitch key) generated independently by each party
# (since it involves their secret keys).
struct CloudKeyPart

    params :: SchemeParameters
    bk_part :: BootstrapKeyPart
    ks :: KeyswitchKey

    function CloudKeyPart(rng, secret_key::SecretKey, shared_key::SharedKey)
        params = secret_key.params
        tlwe_key = TLweKey(rng, tlwe_parameters(params))
        pk = PublicKey(rng, secret_key, tlwe_key, shared_key)
        new(
            params,
            BootstrapKeyPart(rng, secret_key, tlwe_key, shared_key, pk),
            KeyswitchKey(
                rng, params.ks_noise_stddev, keyswitch_parameters(params),
                secret_key.key, tlwe_key))
    end
end


struct MKBootstrapKey
    tgsw_params :: TGswParams
    tlwe_params :: TLweParams
    key :: Array{MKTransformedTGswExpSample, 2}

    function MKBootstrapKey(parts::Array{BootstrapKeyPart, 1})

        parties = length(parts)

        public_keys = [part.public_key for part in parts]

        samples = [
            mk_tgsw_expand(parts[i].key_uni_enc[j], i, public_keys)
            for j in 1:length(parts[1].key_uni_enc), i in 1:parties]

        transformed_samples = forward_transform.(samples)

        new(parts[1].tgsw_params, parts[1].tlwe_params, transformed_samples)
    end
end


struct MKCloudKey

    parties :: Int
    params :: SchemeParameters
    bootstrap_key :: MKBootstrapKey
    keyswitch_key :: Array{KeyswitchKey, 1}

    function MKCloudKey(ck_parts::Array{CloudKeyPart, 1})
        new(
            length(ck_parts),
            ck_parts[1].params,
            MKBootstrapKey([part.bk_part for part in ck_parts]),
            [part.ks for part in ck_parts])
    end
end


function mk_encrypt(rng, secret_keys::Array{SecretKey, 1}, message::Bool)

    mu = encode_message(message ? 1 : -1, 8)

    lwe_params = secret_keys[1].key.params
    alpha = secret_keys[1].params.lwe_noise_stddev
    parties = length(secret_keys)

    a = hcat([rand_uniform_torus32(rng, lwe_params.size) for i in 1:parties]...)
    b = (rand_gaussian_torus32(rng, mu, alpha)
        + reduce(+, a .* hcat([secret_keys[i].key.key for i in 1:parties]...)))

    MKLweSample(lwe_params, a, b, alpha^2)
end


function mk_lwe_phase(sample::MKLweSample, lwe_keys::Array{LweKey, 1})
    parties = length(lwe_keys)
    phases = [
        lwe_phase(LweSample(sample.params, sample.a[:,i], Torus32(0), 0.), lwe_keys[i])
        for i in 1:parties]
    sample.b + reduce(+, phases)
end


function mk_decrypt(secret_keys::Array{SecretKey, 1}, sample::MKLweSample)
    mk_lwe_phase(sample, [sk.key for sk in secret_keys]) > 0
end


mutable struct MKTLweSample
    params :: TLweParams
    a :: Array{TorusPolynomial, 1} # mask (mask_size = 1, so length = parties)
    b :: TorusPolynomial
    current_variance :: Float64

    function MKTLweSample(params::TLweParams, parties::Int)
        p_degree = params.polynomial_degree

        # TODO: assumes mask_size=1 here?
        a = [new_torus_polynomial(p_degree) for i in 1:parties]
        b = new_torus_polynomial(p_degree)

        current_variance = 0.0
        new(params, a, b, current_variance)
    end

    MKTLweSample(params, a, b, cv) = new(params, a, b, cv)
end


Base.:+(x::MKTLweSample, y::MKTLweSample) =
    MKTLweSample(x.params, x.a .+ y.a, x.b + y.b, x.current_variance + y.current_variance)


Base.:-(x::MKTLweSample, y::MKTLweSample) =
    MKTLweSample(x.params, x.a .- y.a, x.b - y.b, x.current_variance + y.current_variance)


# result = (0, ..., 0, mu)
function mk_tlwe_noiseless_trivial(mu::TorusPolynomial, tlwe_params::TLweParams, parties::Int)
    p_degree = tlwe_params.polynomial_degree
    MKTLweSample(
        tlwe_params,
        [torus_polynomial(zeros(Torus32, p_degree)) for i in 1:parties],
        mu,
        0.)
end


function mk_shift_polynomial(sample::MKTLweSample, ai::Int32)
    MKTLweSample(
        sample.params,
        shift_polynomial.(sample.a, ai),
        shift_polynomial(sample.b, ai),
        sample.current_variance)
end


function mk_tgsw_extern_mul(
        sample::MKTLweSample, exp_sample::MKTransformedTGswExpSample, party::Int, parties::Int)

    tlwe_params = sample.params
    tgsw_params = exp_sample.tgsw_params
    p_degree = tlwe_params.polynomial_degree
    decomp_length = tgsw_params.decomp_length

    result = MKTLweSample(tlwe_params, parties)

    dec_a = hcat(decompose.(sample.a, Ref(tgsw_params))...)
    dec_b = decompose(sample.b, tgsw_params)

    tr_dec_a = forward_transform.(dec_a)
    tr_dec_b = forward_transform.(dec_b)

    # c'_i = g^{-1}(a_i)*d1, i<parties, i!=party
    for i in 1:parties
        result.a[i].coeffs .= 0
        if i != party
            for l in 1:decomp_length
                result.a[i] += inverse_transform(tr_dec_a[l,i] * exp_sample.y[l,party])
            end
        end
    end

    # c'_party = \sum g^{-1}(a_i)*yi + g^{-1}(b)*c1
    result.a[party].coeffs .= 0
    for i in 1:parties
        for l in 1:decomp_length
            result.a[party] += inverse_transform(tr_dec_a[l,i] * exp_sample.y[l,i])
        end
    end
    for l in 1:decomp_length
        result.a[party] += inverse_transform(tr_dec_b[l] * exp_sample.c1[l])
    end

    # c'_parties = \sum g^{-1}(a_i)*xi + g^{-1}(b)*c0
    result.b.coeffs .= 0
    for i in 1:parties
        for l in 1:decomp_length
            result.b += inverse_transform(tr_dec_a[l,i] * exp_sample.x[l,i])
        end
    end
    for l in 1:decomp_length
        result.b += inverse_transform(tr_dec_b[l] * exp_sample.c0[l])
    end

    # ORIG_TODO: calculate current_variance

    result
end


function mk_mux_rotate(
        accum::MKTLweSample, sample::MKTransformedTGswExpSample,
        barai::Int32, party::Int, parties::Int)
    # ACC = BKi*[(X^barai-1)*ACC]+ACC
    temp_result = mk_shift_polynomial(accum, barai) - accum
    accum + mk_tgsw_extern_mul(temp_result, sample, party, parties)
end


function mk_blind_rotate(accum::MKTLweSample, bk::MKBootstrapKey, bara::Array{Int32, 2})
    n, parties = size(bara)
    for i in 1:parties
        for j in 1:n
            baraij = bara[j,i]
            if baraij == 0
                continue
            end
            accum = mk_mux_rotate(accum, bk.key[j,i], baraij, i, parties)
        end
    end
    accum
end


function mk_tlwe_extract_sample(x::MKTLweSample)
    # Iterating over parties here, not mask elements! (mask_size = 1)
    # TODO: correct for mask_size > 1
    a = hcat([reverse_polynomial(p).coeffs for p in x.a]...)
    b = x.b.coeffs[1]
    lwe_params = LweParams(x.params.polynomial_degree * x.params.mask_size)
    MKLweSample(lwe_params, a, b, 0.) # TODO: calculate the current variance
end


function mk_blind_rotate_and_extract(
        v::TorusPolynomial, bk::MKBootstrapKey, barb::Int32, bara::Array{Int32, 2})
    parties = size(bara, 2)
    testvectbis = shift_polynomial(v, -barb)
    acc = mk_tlwe_noiseless_trivial(testvectbis, bk.tlwe_params, parties)
    acc = mk_blind_rotate(acc, bk, bara)
    mk_tlwe_extract_sample(acc)
end


function mk_bootstrap_wo_keyswitch(bk::MKBootstrapKey, mu::Torus32, x::MKLweSample)

    p_degree = bk.tlwe_params.polynomial_degree

    barb = decode_message(x.b, p_degree * 2)
    bara = decode_message.(x.a, p_degree * 2)

    #the initial testvec = [mu,mu,mu,...,mu]
    testvect = torus_polynomial(repeat([mu], p_degree))

    mk_blind_rotate_and_extract(testvect, bk, barb, bara)
end


function mk_bootstrap(ck::MKCloudKey, mu::Torus32, x::MKLweSample)
    u = mk_bootstrap_wo_keyswitch(ck.bootstrap_key, mu, x)
    mk_keyswitch(ck.keyswitch_key, u)
end


function mk_lwe_noiseless_trivial(mu::Torus32, params::LweParams, parties::Int)
    MKLweSample(params, zeros(Torus32, params.size, parties), mu, 0.)
end


function mk_gate_nand(ck::MKCloudKey, x::MKLweSample, y::MKLweSample)
    temp = (
        mk_lwe_noiseless_trivial(encode_message(1, 8), x.params, ck.parties)
        - x - y)
    mk_bootstrap(ck, encode_message(1, 8), temp)
end


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
