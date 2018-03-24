push!(LOAD_PATH, "../src")

using TFHE


function int_to_bitarray(x::Int16)
    BitArray([((x >> i) & 1 != 0) for i in 0:15])
end


function bitarray_to_int(x::BitArray)
    int_answer = Int16(0)
    for i in 0:15
        int_answer |= (x[i+1]<<i)
    end
    int_answer
end


function encrypt()

    rng = MersenneTwister(123)
    secret_key, cloud_key = tfhe_key_pair(rng)

    plaintext1 = Int16(2017)
    bits1 = int_to_bitarray(plaintext1)

    plaintext2 = Int16(42)
    bits2 = int_to_bitarray(plaintext2)

    ciphertext1 = tfhe_encrypt(rng, secret_key, bits1)
    ciphertext2 = tfhe_encrypt(rng, secret_key, bits2)

    secret_key, cloud_key, ciphertext1, ciphertext2
end


# elementary full comparator gate that is used to compare the i-th bit:
#   input: ai and bi the i-th bit of a and b
#          lsb_carry: the result of the comparison on the lowest bits
#   algo: if (a==b) return lsb_carry else return b
function encrypted_compare_bit!(cloud_key, result, a, b, lsb_carry, tmp)
    tfhe_gate_XNOR!(cloud_key, tmp, a, b)
    tfhe_gate_MUX!(cloud_key, result, tmp, lsb_carry, a)
end


# this function compares two multibit words, and puts the max in result
function encrypted_minimum!(cloud_key, result, a, b)

    nb_bits = length(result)

    params = tfhe_parameters(cloud_key)

    tmp1 = empty_ciphertext(params, 1)
    tmp2 = empty_ciphertext(params, 1)

    # initialize the carry to 0
    tfhe_gate_CONSTANT!(cloud_key, tmp1, [false])

    # run the elementary comparator gate n times
    for i in 1:length(result)
        encrypted_compare_bit!(cloud_key, tmp1, a[i], b[i], tmp1, tmp2)
    end

    # tmp1 is the result of the comparaison: 0 if a is larger, 1 if b is larger
    # select the max and copy it to the result
    tfhe_gate_MUX!(cloud_key, result, tmp1, b, a)
end


function process(cloud_key, ciphertext1, ciphertext2)

    # if necessary, the params are inside the key
    params = tfhe_parameters(cloud_key)

    # do some operations on the ciphertexts: here, we will compute the
    # minimum of the two
    result = empty_ciphertext(params, size(ciphertext1)...)
    encrypted_minimum!(cloud_key, result, ciphertext1, ciphertext2)

    result
end


function verify(secret_key, answer)
    answer_bits = tfhe_decrypt(secret_key, answer)
    int_answer = bitarray_to_int(answer_bits)
    println("Answer: $int_answer")
end


secret_key, cloud_key, ciphertext1, ciphertext2 = encrypt()
answer = process(cloud_key, ciphertext1, ciphertext2)
verify(secret_key, answer)
