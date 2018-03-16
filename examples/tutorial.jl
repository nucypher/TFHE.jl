push!(LOAD_PATH, "../src")

using TFHE


function encrypt()
    # generate a keyset
    minimum_lambda = 110
    params = new_default_gate_bootstrapping_parameters(minimum_lambda)

    # generate a random key
    seed = UInt32[314, 1592, 657]
    tfhe_random_generator_setSeed(seed)

    key = new_random_gate_bootstrapping_secret_keyset(params)

    # generate encrypt the 16 bits of 2017
    plaintext1 = Int16(2017)
    ciphertext1 = new_gate_bootstrapping_ciphertext_array(16, params)
    for i in 0:15
        bootsSymEncrypt(ciphertext1[i+1], ((plaintext1>>i) & 1) != 0, key)
    end

    # generate encrypt the 16 bits of 42
    plaintext2 = Int16(42)
    ciphertext2 = new_gate_bootstrapping_ciphertext_array(16, params)
    for i in 0:15
        bootsSymEncrypt(ciphertext2[i+1], ((plaintext2>>i) & 1) != 0, key)
    end


    key, key.cloud, ciphertext1, ciphertext2
end


# elementary full comparator gate that is used to compare the i-th bit:
#   input: ai and bi the i-th bit of a and b
#          lsb_carry: the result of the comparison on the lowest bits
#   algo: if (a==b) return lsb_carry else return b
function compare_bit(
        result::LweSample, a::LweSample, b::LweSample, lsb_carry::LweSample,
        tmp::LweSample, bk::TFheGateBootstrappingCloudKeySet)
    bootsXNOR(tmp, a, b, bk)
    bootsMUX(result, tmp, lsb_carry, a, bk)
end

# this function compares two multibit words, and puts the max in result
function minimum(
        result::Array{LweSample, 1}, a::Array{LweSample, 1}, b::Array{LweSample, 1}, nb_bits::Int,
        bk::TFheGateBootstrappingCloudKeySet)

    tmps = new_gate_bootstrapping_ciphertext_array(2, bk.params)

    # initialize the carry to 0
    bootsCONSTANT(tmps[1], Int32(0), bk)
    # run the elementary comparator gate n times
    for i in 0:(nb_bits-1)
        compare_bit(tmps[1], a[i+1], b[i+1], tmps[0+1], tmps[1+1], bk)
    end

    # tmps[0] is the result of the comparaison: 0 if a is larger, 1 if b is larger
    # select the max and copy it to the result
    for i in 0:(nb_bits-1)
        bootsMUX(result[i+1], tmps[0+1], b[i+1], a[i+1], bk)
    end
end


function process(cloud_key, ciphertext1, ciphertext2)

    # if necessary, the params are inside the key
    params = cloud_key.params

    # do some operations on the ciphertexts: here, we will compute the
    # minimum of the two
    result = new_gate_bootstrapping_ciphertext_array(16, params)
    minimum(result, ciphertext1, ciphertext2, 16, cloud_key)

    result
end


function verify(secret_key, answer)

    # if necessary, the params are inside the key
    params = secret_key.params

    # decrypt and rebuild the 16-bit plaintext answer
    int_answer = Int16(0)
    for i in 0:15
        ai = bootsSymDecrypt(answer[i+1], secret_key)
        int_answer |= (ai<<i)
    end

    println("Answer: $int_answer")
end


secret_key, cloud_key, ciphertext1, ciphertext2 = encrypt()
answer = process(cloud_key, ciphertext1, ciphertext2)
verify(secret_key, answer)
