# LWE

"""
A structure representing an encrypted bit in multi-key TFHE.
"""
struct MKLweSample

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


function mk_lwe_phase(sample::MKLweSample, lwe_keys::Array{LweKey, 1})
    parties = length(lwe_keys)
    phases = [
        lwe_phase(LweSample(sample.params, sample.a[:,i], Torus32(0), 0.), lwe_keys[i])
        for i in 1:parties]
    sample.b + reduce(+, phases)
end


function mk_lwe_noiseless_trivial(mu::Torus32, params::LweParams, parties::Int)
    MKLweSample(params, zeros(Torus32, params.size, parties), mu, 0.)
end


# TLWE


struct MKTLweSample
    params :: TLweParams
    a :: Array{TorusPolynomial, 1} # mask (mask_size = 1, so length = parties)
    b :: TorusPolynomial
    current_variance :: Float64

    function MKTLweSample(
            params::TLweParams, a::Array{<:TorusPolynomial, 1}, b::TorusPolynomial,
            current_variance::Float64)
        new(params, a, b, current_variance)
    end
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


function mk_mul_by_monomial(sample::MKTLweSample, ai::Int32)
    MKTLweSample(
        sample.params,
        mul_by_monomial.(sample.a, ai),
        mul_by_monomial(sample.b, ai),
        sample.current_variance)
end


function mk_tlwe_extract_sample(x::MKTLweSample)
    # Iterating over parties here, not mask elements! (mask_size = 1)
    # TODO: (issue #2) correct for mask_size > 1
    a = hcat([reverse_polynomial(p).coeffs for p in x.a]...)
    b = x.b.coeffs[1]
    lwe_params = LweParams(x.params.polynomial_degree * x.params.mask_size)
    MKLweSample(lwe_params, a, b, 0.) # TODO: (issue #7) calculate the current variance
end


# Internal keys


struct SharedKey
    tgsw_params :: TGswParams
    tlwe_params :: TLweParams
    a :: Array{TorusPolynomial, 1}

    function SharedKey(rng, tgsw_params::TGswParams, tlwe_params::TLweParams)
        decomp_length = tgsw_params.decomp_length
        p_degree = tlwe_params.polynomial_degree
        a = [torus_polynomial(rand_uniform_torus32(rng, p_degree)) for i in 1:decomp_length]
        new(tgsw_params, tlwe_params, a)
    end
end


struct PublicKey
    tgsw_params :: TGswParams
    tlwe_params :: TLweParams
    b :: Array{TorusPolynomial, 1}

    function PublicKey(
            rng, tlwe_key::TLweKey, alpha::Float64, shared::SharedKey, tgsw_params::TGswParams)

        tlwe_params = tlwe_key.params
        decomp_length = tgsw_params.decomp_length
        p_degree = tlwe_params.polynomial_degree

        # The right part of the public key is `b_i = e_i + a*s_i`,
        # where `a` is shared between the parties
        # TODO: (issue #2) [1] was omitted in the original!
        #       It works while k=1, but will fail otherwise
        # TODO: (issue #2) this is basically tgsw_encrypt_zero() for mask_size=1
        b = [(
                transformed_mul(tlwe_key.key[1], shared.a[i])
                + torus_polynomial(rand_gaussian_torus32(rng, zero(Int32), alpha, p_degree)))
            for i in 1:decomp_length]

        new(tgsw_params, tlwe_params, b)
    end
end


# TGSW


# Uni-encrypted TGSW sample
struct MKTGswUESample

    tgsw_params :: TGswParams
    tlwe_params :: TLweParams

    c0 :: Array{TorusPolynomial, 1}
    c1 :: Array{TorusPolynomial, 1}
    d0 :: Array{TorusPolynomial, 1}
    d1 :: Array{TorusPolynomial, 1}
    f0 :: Array{TorusPolynomial, 1}
    f1 :: Array{TorusPolynomial, 1}

    current_variance :: Float64 # avg variance of the sample

    function MKTGswUESample(
            tgsw_params::TGswParams, tlwe_params::TLweParams,
            c0::Array{<: TorusPolynomial, 1}, c1::Array{<: TorusPolynomial, 1},
            d0::Array{<: TorusPolynomial, 1}, d1::Array{<: TorusPolynomial, 1},
            f0::Array{<: TorusPolynomial, 1}, f1::Array{<: TorusPolynomial, 1},
            current_variance::Float64)

        decomp_length = tgsw_params.decomp_length
        p_degree = tlwe_params.polynomial_degree

        @assert size(c0) == (decomp_length,)
        @assert size(c1) == (decomp_length,)
        @assert size(d0) == (decomp_length,)
        @assert size(d1) == (decomp_length,)
        @assert size(f0) == (decomp_length,)
        @assert size(f1) == (decomp_length,)

        new(tgsw_params, tlwe_params, c0, c1, d0, d1, f0, f1, current_variance)
    end
end


# Encrypt an integer value
# In the paper: RGSW.UniEnc
# Similar to tgsw_encrypt()/tlwe_encrypt(), except the public key is supplied externally.
function mk_tgsw_encrypt(
        rng, message::Int32, alpha::Float64,
        tlwe_key::TLweKey, shared_key::SharedKey, public_key::PublicKey)

    tgsw_params = shared_key.tgsw_params
    tlwe_params = tlwe_key.params
    p_degree = tlwe_params.polynomial_degree
    decomp_length = tgsw_params.decomp_length

    # The shared randomness
    r = int_polynomial(rand_uniform_bool(rng, p_degree))

    # C = (c0,c1) \in T^2dg, with c0 = s_party*c1 + e_c + m*g
    c1 = [torus_polynomial(rand_uniform_torus32(rng, p_degree)) for i in 1:decomp_length]
    # TODO: (issue #2) it was just key, not key[1] in the original. Hardcoded mask_size=1? Check!
    c0 = [
            (torus_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, p_degree))
            + transformed_mul(tlwe_key.key[1], c1[i])
            + message * tgsw_params.gadget_values[i])
            for i in 1:decomp_length]

    # D = (d0, d1) = r*[Pkey_party | Pkey_parties] + [E0|E1] + [0|m*g] \in T^2dg
    d1 = [
            (torus_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, p_degree))
            + transformed_mul(r, shared_key.a[i])
            + message * tgsw_params.gadget_values[i])
            for i in 1:decomp_length]
    d0 = [
            (torus_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, p_degree))
            + transformed_mul(r, public_key.b[i]))
            for i in 1:decomp_length]

    # F = (f0,f1) \in T^2dg, with f0 = s_party*f1 + e_f + r*g
    f1 = [torus_polynomial(rand_uniform_torus32(rng, p_degree)) for i in 1:decomp_length]
    # TODO: (issue #2) it was just key, not key[1] in the original. Hardcoded mask_size=1? Check!
    f0 = [
            (torus_polynomial(rand_gaussian_torus32(rng, Int32(0), alpha, p_degree))
            + transformed_mul(tlwe_key.key[1], f1[i])
            + r * tgsw_params.gadget_values[i])
            for i in 1:decomp_length]

    MKTGswUESample(tgsw_params, tlwe_params, c0, c1, d0, d1, f0, f1, alpha^2)
end


# TGSW sample expanded for all the parties.
# The full matrix has size `(parties+1)x(parties+1)` and its elements are
# `decomp_length`-vectors of `TorusPolynomial`.
# Since the C matrix returned by RGSW.Expand in the paper is sparse and has repeating elements,
# we are keeping only nonzero elements of it, namely:
# To save memory, we are only keeping the distinct elements of it, namely:
# c0 == C_{1,1},
# c1 == C_{1,party+1}
# x == C_{i,1}, i = 2...parties+1
# y == C_{i,party+1}, i == 2...parties+1
struct MKTGswExpSample

    tgsw_params :: TGswParams
    tlwe_params :: TLweParams

    x :: Array{TorusPolynomial, 2}
    y :: Array{TorusPolynomial, 2}
    c0 :: Array{TorusPolynomial, 1}
    c1 :: Array{TorusPolynomial, 1}

    current_variance :: Float64 # avg variance of the sample

    function MKTGswExpSample(
            tgsw_params::TGswParams, tlwe_params::TLweParams,
            x::Array{<: TorusPolynomial, 2},
            y::Array{<: TorusPolynomial, 2},
            c0::Array{<: TorusPolynomial, 1},
            c1::Array{<: TorusPolynomial, 1},
            current_variance::Float64)

        decomp_length = tgsw_params.decomp_length
        p_degree = tlwe_params.polynomial_degree
        parties = size(x, 2)

        @assert size(x) == (decomp_length, parties)
        @assert size(y) == (decomp_length, parties)
        @assert size(c0) == (decomp_length,)
        @assert size(c1) == (decomp_length,)

        new(tgsw_params, tlwe_params, x, y, c0, c1, current_variance)
    end
end


struct MKTransformedTGswExpSample

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


# In the paper: RGSW.Expand
function mk_tgsw_expand(sample::MKTGswUESample, party::Int, public_keys::Array{PublicKey, 1})

    tlwe_params = sample.tlwe_params
    tgsw_params = sample.tgsw_params
    p_degree = sample.tlwe_params.polynomial_degree
    decomp_length = tgsw_params.decomp_length
    parties = length(public_keys)

    # Calculate the distinct elements of the expanded matrix C.
    # See the description of `MKTGswUESample` for details.

    zero_p = torus_polynomial(zeros(Torus32, p_degree))

    # g^{-1}(b_i[j] - b_party[j]) = [u_0, ...,u_dg-1]
    decomp = [
        (i == party ?
            nothing : # unused
            decompose(public_keys[i].b[j] - public_keys[party].b[j], tgsw_params)
            )
        for j in 1:decomp_length, i in 1:parties
        ]

    x = [
        (sample.d0[j] + (i == party ?
            zero_p :
            # add xi[j] = <g^{-1}(b_i[j] - b_party[j]), f0>
            sum(transformed_mul(decomp[j,i][l], sample.f0[l]) for l in 1:decomp_length)))
        for j in 1:decomp_length, i in 1:parties
        ]

    y = [
        (i == party ?
            sample.d1[j] :
            # add yi[j] = <g^{-1}(b_i[j] - b_party[j]), f1>
            sum(transformed_mul(decomp[j,i][l], sample.f1[l]) for l in 1:decomp_length)
            )
        for j in 1:decomp_length, i in 1:parties
        ]

    # TODO: (issue #7) calculate the current variance correctly
    MKTGswExpSample(tgsw_params, tlwe_params, x, y, sample.c0, sample.c1, 0.0)
end


function mk_tgsw_extern_mul(
        sample::MKTLweSample, exp_sample::MKTransformedTGswExpSample, party::Int, parties::Int)

    tlwe_params = sample.params
    tgsw_params = exp_sample.tgsw_params
    p_degree = tlwe_params.polynomial_degree
    decomp_length = tgsw_params.decomp_length

    dec_a = hcat(decompose.(sample.a, Ref(tgsw_params))...)
    dec_b = decompose(sample.b, tgsw_params)

    # It would be faster to do all the summations in transformed space and then transform back.
    # But, because of the transform used (FFT), the precision is limited,
    # and we are close to that limit.
    # Elements in `dec_a` and `dec_b` use `bs_log2_base` bits.
    # Elements in `exp_sample` arrays have full range of Int32 (32 bits).
    # So if we want to do all the calculations in this function in transformed space, that takes
    # `bs_log2_base + 32 + log2(p_degree) + ceil(log2(decomp_length * (parties+1)))` bits
    # == 8 + 32 + 10 + ceil(log2(4 * 3)) = 54 bits - not good.

    tr_dec_a = forward_transform.(dec_a)
    tr_dec_b = forward_transform.(dec_b)

    # c'_i = g^{-1}(a_i)*d1, i<parties, i!=party
    # c'_party = \sum g^{-1}(a_i)*yi + g^{-1}(b)*c1
    a = [
        (i == party ?
            (sum([inverse_transform(tr_dec_a[l,j] * exp_sample.y[l,j])
                for l in 1:decomp_length, j in 1:parties])
            + sum([inverse_transform(tr_dec_b[l] * exp_sample.c1[l])
                for l in 1:decomp_length])) :
            sum([inverse_transform(tr_dec_a[l,i] * exp_sample.y[l,party])
                for l in 1:decomp_length]))
        for i in 1:parties]

    # c'_parties = \sum g^{-1}(a_i)*xi + g^{-1}(b)*c0
    b = (sum([inverse_transform(tr_dec_a[l,i] * exp_sample.x[l,i])
            for l in 1:decomp_length, i in 1:parties])
        + sum([inverse_transform(tr_dec_b[l] * exp_sample.c0[l])
            for l in 1:decomp_length]))

    # TODO: (issue #7) calculate current_variance
    MKTLweSample(tlwe_params, a, b, 0.)
end


# Keyswitch


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


# Bootstrap


# A part of the bootstrap key generated independently by each party
# (since it involves their secret keys).
struct BootstrapKeyPart

    tgsw_params :: TGswParams
    tlwe_params :: TLweParams
    key_uni_enc :: Array{MKTGswUESample, 1}
    public_key :: PublicKey

    function BootstrapKeyPart(
            rng, lwe_key::LweKey, tlwe_key::TLweKey, alpha::Float64,
            shared_key::SharedKey, public_key::PublicKey)

        tgsw_params = shared_key.tgsw_params
        tlwe_params = shared_key.tlwe_params
        in_out_params = lwe_key.params
        uni_enc = [
            mk_tgsw_encrypt(rng, lwe_key.key[j], alpha, tlwe_key, shared_key, public_key)
            for j in 1:in_out_params.size]

        new(tgsw_params, tlwe_params, uni_enc, public_key)
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


function mk_mux_rotate(
        accum::MKTLweSample, sample::MKTransformedTGswExpSample,
        barai::Int32, party::Int, parties::Int)
    # ACC = BKi*[(X^barai-1)*ACC]+ACC
    temp_result = mk_mul_by_monomial(accum, barai) - accum
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


function mk_blind_rotate_and_extract(
        v::TorusPolynomial, bk::MKBootstrapKey, barb::Int32, bara::Array{Int32, 2})
    parties = size(bara, 2)
    testvectbis = mul_by_monomial(v, -barb)
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


function mk_bootstrap(bk::MKBootstrapKey, ks::Array{KeyswitchKey, 1}, mu::Torus32, x::MKLweSample)
    u = mk_bootstrap_wo_keyswitch(bk, mu, x)
    mk_keyswitch(ks, u)
end
