#=
Homomorphic gates.
Take LWE samples with message in {-1/8, 1/8}, noise<1/16,
return a (bootstrapped) LWE sample with message in {-1/8, 1/8}, noise<1/16
(that is, a negative phase encodes `false`, a positive phase encodes `true`).
=#


function gate_nand(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(1, 8), in_out_params) - x - y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


function gate_or(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(1, 8), in_out_params) + x + y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


function gate_and(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(-1, 8), in_out_params) + x + y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


function gate_xor(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(1, 4), in_out_params) + (x + y) * 2
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


function gate_xnor(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(-1, 4), in_out_params) - (x + y) * 2
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


function gate_not(ck::CloudKey, x::LweSample)
    # Not bootstrapped, the cloud key is just for the sake of interface uniformity.
    -x
end


function gate_constant(ck::CloudKey, value::Bool)
    in_out_params = ck.params.in_out_params
    lwe_noiseless_trivial(encode_message(value ? 1 : -1, 8), in_out_params)
end


function gate_nor(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(-1, 8), in_out_params) - x - y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
ANDNY(x, y) == AND(NOT(x), y)
"""
function gate_andny(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(-1, 8), in_out_params) - x + y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
ANDYN(x, y) == AND(x, NOT(y))
"""
function gate_andyn(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(-1, 8), in_out_params) + x - y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
ORNY(x, y) == OR(NOT(x), y)
"""
function gate_orny(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(1, 8), in_out_params) - x + y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
ORYN(x, y) == OR(x, NOT(y))
"""
function gate_oryn(ck::CloudKey, x::LweSample, y::LweSample)
    in_out_params = ck.params.in_out_params
    result = lwe_noiseless_trivial(encode_message(1, 8), in_out_params) + x - y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
MUX(x, y, z) == x ? y : z == OR(AND(x, y), AND(NOT(x), z))
"""
function gate_mux(ck::CloudKey, x::LweSample, y::LweSample, z::LweSample)

    in_out_params = ck.params.in_out_params
    extracted_params = ck.params.tgsw_params.tlwe_params.extracted_lweparams

    # compute `AND(x, y)`
    t1 = lwe_noiseless_trivial(encode_message(-1, 8), in_out_params) + x + y
    u1 = bootstrap_wo_keyswitch(ck.bootstrap_key, encode_message(1, 8), t1)

    # compute `AND(NOT(x), z)`
    t2 = lwe_noiseless_trivial(encode_message(-1, 8), in_out_params) - x + z
    u2 = bootstrap_wo_keyswitch(ck.bootstrap_key, encode_message(1, 8), t2)

    # compute `OR(u1,u2)`
    t3 = lwe_noiseless_trivial(encode_message(1, 8), extracted_params) + u1 + u2

    keyswitch(ck.keyswitch_key, t3)
end
