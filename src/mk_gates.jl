function mk_gate_nand(ck::MKCloudKey, x::MKLweSample, y::MKLweSample)
    temp = (
        mk_lwe_noiseless_trivial(encode_message(1, 8), x.params, ck.parties)
        - x - y)
    mk_bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), temp)
end
