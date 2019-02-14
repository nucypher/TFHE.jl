#*#*****************************************
# zones on the torus -> to see
#*#*****************************************


#=
 * Homomorphic bootstrapped NAND gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_NAND(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    #compute: (0,1/8) - ca - cb
    NandConst = modSwitchToTorus32(1, 8)
    temp_result = lweNoiselessTrivial(NandConst, in_out_params)
    temp_result -= ca
    temp_result -= cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped OR gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_OR(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    #compute: (0,1/8) + ca + cb
    OrConst = modSwitchToTorus32(1, 8)
    temp_result = lweNoiselessTrivial(OrConst, in_out_params)
    temp_result += ca
    temp_result += cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped AND gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_AND(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    #compute: (0,-1/8) + ca + cb
    AndConst = modSwitchToTorus32(-1, 8)
    temp_result = lweNoiselessTrivial(AndConst, in_out_params)
    temp_result += ca
    temp_result += cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped XOR gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_XOR(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    #compute: (0,1/4) + 2*(ca + cb)
    XorConst = modSwitchToTorus32(1, 4)
    temp_result = lweNoiselessTrivial(XorConst, in_out_params)
    temp_result += ca * 2
    temp_result += cb * 2

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped XNOR gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_XNOR(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    result = TFHEEncryptedBit(in_out_params)
    temp_result = TFHEEncryptedBit(in_out_params)

    #compute: (0,-1/4) + 2*(-ca-cb)
    XnorConst = modSwitchToTorus32(-1, 4)
    temp_result = lweNoiselessTrivial(XnorConst, in_out_params)
    temp_result -= ca * 2
    temp_result -= cb * 2

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped NOT gate (doesn't need to be bootstrapped)
 * Takes in input 1 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_NOT(bk::TFHECloudKey, ca::TFHEEncryptedBit)
    -ca
end


#=
 * Homomorphic Trivial Constant gate (doesn't need to be bootstrapped)
 * Takes a boolean value)
 * Outputs a LWE sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_CONSTANT(bk::TFHECloudKey, value::Bool)
    in_out_params = bk.params.in_out_params
    result = TFHEEncryptedBit(in_out_params)
    MU = modSwitchToTorus32(1, 8)
    lweNoiselessTrivial(value ? MU : -MU, in_out_params)
end


#=
 * Homomorphic bootstrapped NOR gate
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_NOR(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    #compute: (0,-1/8) - ca - cb
    NorConst = modSwitchToTorus32(-1, 8)
    temp_result = lweNoiselessTrivial(NorConst, in_out_params)
    temp_result -= ca
    temp_result -= cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped AndNY Gate: not(a) and b
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_ANDNY(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    #compute: (0,-1/8) - ca + cb
    AndNYConst = modSwitchToTorus32(-1, 8)
    temp_result = lweNoiselessTrivial(AndNYConst, in_out_params)
    temp_result -= ca
    temp_result += cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped AndYN Gate: a and not(b)
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_ANDYN(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    #compute: (0,-1/8) + ca - cb
    AndYNConst = modSwitchToTorus32(-1, 8)
    temp_result = lweNoiselessTrivial(AndYNConst, in_out_params)
    temp_result += ca
    temp_result -= cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped OrNY Gate: not(a) or b
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_ORNY(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    #compute: (0,1/8) - ca + cb
    OrNYConst = modSwitchToTorus32(1, 8)
    temp_result = lweNoiselessTrivial(OrNYConst, in_out_params)
    temp_result -= ca
    temp_result += cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped OrYN Gate: a or not(b)
 * Takes in input 2 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_ORYN(
        bk::TFHECloudKey, ca::TFHEEncryptedBit, cb::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params

    #compute: (0,1/8) + ca - cb
    OrYNConst = modSwitchToTorus32(1, 8)
    temp_result = lweNoiselessTrivial(OrYNConst, in_out_params)
    temp_result += ca
    temp_result -= cb

    #if the phase is positive, the result is 1/8
    #if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(bk.bkFFT, MU, temp_result)
end


#=
 * Homomorphic bootstrapped Mux(a,b,c) = a?b:c = a*b + not(a)*c
 * Takes in input 3 LWE samples (with message space [-1/8,1/8], noise<1/16)
 * Outputs a LWE bootstrapped sample (with message space [-1/8,1/8], noise<1/16)
=#
function tfhe_gate_MUX(
        bk::TFHECloudKey,
        a::TFHEEncryptedBit, b::TFHEEncryptedBit, c::TFHEEncryptedBit)

    MU = modSwitchToTorus32(1, 8)
    in_out_params = bk.params.in_out_params
    extracted_params = bk.params.tgsw_params.tlwe_params.extracted_lweparams

    #compute "AND(a,b)": (0,-1/8) + a + b
    AndConst = modSwitchToTorus32(-1, 8)
    temp_result = lweNoiselessTrivial(AndConst, in_out_params)
    temp_result += a
    temp_result += b
    # Bootstrap without KeySwitch
    u1 = tfhe_bootstrap_woKS_FFT(bk.bkFFT, MU, temp_result)


    #compute "AND(not(a),c)": (0,-1/8) - a + c
    temp_result = lweNoiselessTrivial(AndConst, in_out_params)
    temp_result -= a
    temp_result += c
    # Bootstrap without KeySwitch
    u2 = tfhe_bootstrap_woKS_FFT(bk.bkFFT, MU, temp_result)

    # Add u1=u1+u2
    MuxConst = modSwitchToTorus32(1, 8)
    temp_result1 = lweNoiselessTrivial(MuxConst, extracted_params)
    temp_result1 += u1
    temp_result1 += u2
    # Key switching
    lweKeySwitch(bk.bkFFT.ks, temp_result1)
end
