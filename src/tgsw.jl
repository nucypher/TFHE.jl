struct TGswParams
    decomp_length :: Int32 # decomposition length
    log2_base :: Int32 # log2(decomposition base)
    tlwe_params :: TLweParams # Params of each row

    gadget_values :: Array{Torus32, 1} # powers of the decomposition base
    offset :: Int32 # decomposition offset

    function TGswParams(decomp_length::Int, log2_base::Int, tlwe_params::TLweParams)

        decomp_range = 1:decomp_length

        # Nonzero entries of the gadget matrix.
        # 1/(base^(i+1)) as a Torus32 (the values of the gadget matrix)
        gadget_values = Torus32(2) .^ (32 .- decomp_range .* log2_base)

        # A precalculated offset value for decomposition of torus elements.
        # offset = base/2 * Sum{j=1..decomp_length} 2^(32 - j * bs_log2_base)
        offset = signed(UInt32(sum(gadget_values) * (2^(log2_base - 1))))

        new(decomp_length, log2_base, tlwe_params, gadget_values, offset)
    end
end


struct TGswKey
    params :: TGswParams
    tlwe_key :: TLweKey

    function TGswKey(rng::AbstractRNG, params::TGswParams)
        tlwe_key = TLweKey(rng, params.tlwe_params)
        new(params, tlwe_key)
    end
end


struct TGswSample
    samples :: Array{TLweSample, 2}

    function TGswSample(params::TGswParams)
        mask_size = params.tlwe_params.mask_size
        l = params.decomp_length
        samples = [TLweSample(params.tlwe_params) for i in 1:((mask_size + 1) * l)]
        new(reshape(samples, Int64(l), mask_size + 1))
    end

    TGswSample(samples) = new(samples)
end


struct TransformedTGswSample
    samples:: Array{TransformedTLweSample, 2}

    function TransformedTGswSample(params::TGswParams)
        mask_size = params.tlwe_params.mask_size
        l = params.decomp_length
        samples = [TransformedTLweSample(params.tlwe_params) for i in 1:((mask_size + 1) * l)]
        new(reshape(samples, Int64(l), mask_size + 1))
    end

    TransformedTGswSample(samples) = new(samples)
end


"""
Returns `sample + message * h`, where `h` is the gadget matrix
(`h_ijk = delta_jk / B^i`, `i=1...decomp_length`, `j,k=1...mask_size+1`,
`B` is the decomposition base).
The dimensions of the TLWE sample array in the TGSW sample correspond to the first two indices
(`i` and `j`), and the index `k` corresponds to the vector of polynomials in a single TLWE sample.
"""
function tgsw_add_gadget_times_message(sample::TGswSample, message::Int32, params::TGswParams)
    mask_size = params.tlwe_params.mask_size
    decomp_length = params.decomp_length
    gadget = params.gadget_values

    # Since the gadget matrix is block-diagonal, we avoid extra work by only doing the addition
    # where its elements are nonzero.
    # (A violation of the Liskov principle here, which could be avoided by using special
    # types for zero TLWE samples and zero polynomials, but that's just too much infrastructure
    # for one small function)
    TGswSample([
        TLweSample([
            j == k ? sample.samples[i,j].a[k] + message * gadget[i] : sample.samples[i,j].a[k]
            for k in 1:mask_size+1],
            sample.samples[i,j].current_variance) # TODO: calculate current_variance correctly
        for i in 1:decomp_length, j in 1:mask_size+1])

end


function tgsw_encrypt_zero(rng::AbstractRNG, alpha::Float64, key::TGswKey)
    params = key.params
    mask_size = params.tlwe_params.mask_size
    decomp_length = params.decomp_length
    TGswSample([
        tlwe_encrypt_zero(rng, alpha, key.tlwe_key, params.tlwe_params)
        for i in 1:decomp_length, j in 1:mask_size+1])
end


function tgsw_encrypt(rng::AbstractRNG, message::Int32, alpha::Float64, key::TGswKey)
    tgsw_zero = tgsw_encrypt_zero(rng, alpha, key)
    tgsw_add_gadget_times_message(tgsw_zero, message, key.params)
end


"""
Given the decomposition length `l`, decompose each coefficient of the given polynomial
into `l` values and store them in `l` polynomials.

For each value `x` in the real torus, the decomposition procedure floors it
to a multiple of `1/B^l` (where `B == 2^log2_base` is the decomposition base)
and finds `l` values `a_i` in `[-B/2, B/2) such that `x = sum(a_i / B^i, i=1...l)`.
"""
function decompose(sample::TorusPolynomial, params::TGswParams)

    decomp_length = params.decomp_length
    log2_base = params.log2_base

    mask = Int32((1 << params.log2_base) - 1)
    part_offset = Int32(1 << (params.log2_base - 1))
    offset = params.offset

    # Since we want results in the range `[-B/2, B/2)` instead of `[0, B)`,
    # we add an `offset = sum(base powers)`,
    # decompose normally by calculating remainders
    # (since our base powers are powers of 2, these are implemented as shifts and ANDs)
    # and subtract back part of the offset (`B/2`).

    [int_polynomial(
        @. (((sample.coeffs + offset) >> (32 - power * log2_base)) & mask) - part_offset)
        for power in 1:decomp_length]
end


forward_transform(source::TGswSample) =
    TransformedTGswSample(forward_transform.(source.samples))


# External product (*): accum = gsw (*) accum
function tgsw_extern_mul(accum::TLweSample, gsw::TransformedTGswSample, params::TGswParams)
    dec_accum = hcat(decompose.(accum.a, Ref(params))...)
    tr_dec_accum = forward_transform.(dec_accum)
    inverse_transform(sum(gsw.samples .* tr_dec_accum))
end
