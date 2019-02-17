#*#*****************************************
# zones on the torus -> to see
#*#*****************************************


#=
 * Homomorphic bootstrapped NAND gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_NAND(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    #compute: (0,1/8) - ca - cb
    NandConst = modSwitchToTorus32(1, 8)
    temp_result = lwe_noiseless_trivial(NandConst, in_out_params)
    temp_result -= ca
    temp_result -= cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped OR gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_OR(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    #compute: (0,1/8) + ca + cb
    OrConst = modSwitchToTorus32(1, 8)
    temp_result = lwe_noiseless_trivial(OrConst, in_out_params)
    temp_result += ca
    temp_result += cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped AND gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_AND(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    #compute: (0,-1/8) + ca + cb
    AndConst = modSwitchToTorus32(-1, 8)
    temp_result = lwe_noiseless_trivial(AndConst, in_out_params)
    temp_result += ca
    temp_result += cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped XOR gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_XOR(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    #compute: (0,1/4) + 2*(ca + cb)
    XorConst = modSwitchToTorus32(1, 4)
    temp_result = lwe_noiseless_trivial(XorConst, in_out_params)
    temp_result += ca * 2
    temp_result += cb * 2

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped XNOR gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_XNOR(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    result = TFHEEncryptedBit(in_out_params)
    temp_result = TFHEEncryptedBit(in_out_params)

    #compute: (0,-1/4) + 2*(-ca-cb)
    XnorConst = modSwitchToTorus32(-1, 4)
    temp_result = lwe_noiseless_trivial(XnorConst, in_out_params)
    temp_result -= ca * 2
    temp_result -= cb * 2

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped NOT gate (doesn't need to be bootstrapped)
 * Takes in input 1 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_NOT(ck::TFHECloudKey, ca::TFHEEncryptedBit)
    -ca
end


#=
 * Homomorphic Trivial Constant gate (doesn't need to be bootstrapped)
 * Takes a boolean value)
 * Outputs a LWE sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_CONSTANT(ck::TFHECloudKey, value::Bool)
    in_out_params = ck.params.in_out_params
    result = TFHEEncryptedBit(in_out_params)
    MU = modSwitchToTorus32(1, 8)
    lwe_noiseless_trivial(value ? MU : -MU, in_out_params)
end


#=
 * Homomorphic bootstrapped NOR gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_NOR(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    #compute: (0,-1/8) - ca - cb
    NorConst = modSwitchToTorus32(-1, 8)
    temp_result = lwe_noiseless_trivial(NorConst, in_out_params)
    temp_result -= ca
    temp_result -= cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped AndNY Gate: not(a) and b
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_ANDNY(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    #compute: (0,-1/8) - ca + cb
    AndNYConst = modSwitchToTorus32(-1, 8)
    temp_result = lwe_noiseless_trivial(AndNYConst, in_out_params)
    temp_result -= ca
    temp_result += cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped AndYN Gate: a and not(b)
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_ANDYN(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    #compute: (0,-1/8) + ca - cb
    AndYNConst = modSwitchToTorus32(-1, 8)
    temp_result = lwe_noiseless_trivial(AndYNConst, in_out_params)
    temp_result += ca
    temp_result -= cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped OrNY Gate: not(a) or b
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_ORNY(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    #compute: (0,1/8) - ca + cb
    OrNYConst = modSwitchToTorus32(1, 8)
    temp_result = lwe_noiseless_trivial(OrNYConst, in_out_params)
    temp_result -= ca
    temp_result += cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped OrYN Gate: a or not(b)
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_ORYN(
        ck::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params

    #compute: (0,1/8) + ca - cb
    OrYNConst = modSwitchToTorus32(1, 8)
    temp_result = lwe_noiseless_trivial(OrYNConst, in_out_params)
    temp_result += ca
    temp_result -= cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(ck.bootstrap_key, ck.keyswitch_key, MU, temp_result)
end


#=
 * Homomorphic bootstrapped Mux(a,b,c) = a?b:c = a*b + not(a)*c
 * Takes in input 3 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_MUX(
        ck::TFHECloudKey,
        a::TFHEEncryptedBit, b::TFHEEncryptedBit, c::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = ck.params.in_out_params
    extracted_params = ck.params.tgsw_params.tlwe_params.extracted_lweparams

    #compute "AND(a,b)": (0,-1/8) + a + b
    AndConst = modSwitchToTorus32(-1, 8)
    temp_result = lwe_noiseless_trivial(AndConst, in_out_params)
    temp_result += a
    temp_result += b
    # Bootstrap without KeySwitch
    u1 = tfhe_bootstrap_woKS_FFT(ck.bootstrap_key, MU, temp_result)


    #compute "AND(not(a),c)": (0,-1/8) - a + c
    temp_result = lwe_noiseless_trivial(AndConst, in_out_params)
    temp_result -= a
    temp_result += c
    # Bootstrap without KeySwitch
    u2 = tfhe_bootstrap_woKS_FFT(ck.bootstrap_key, MU, temp_result)

    # Add u1=u1+u2
    MuxConst = modSwitchToTorus32(1, 8)
    temp_result1 = lwe_noiseless_trivial(MuxConst, extracted_params)
    temp_result1 += u1
    temp_result1 += u2
    # Key switching
    lweKeySwitch(ck.keyswitch_key, temp_result1)
end
