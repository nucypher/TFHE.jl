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

#=
# elementary full comparator gate that is used to compare the i-th bit:
#   input: ai and bi the i-th bit of a and b
#          lsb_carry: the result of the comparison on the lowest bits
#   algo: if (a==b) return lsb_carry else return b
function encrypted_compare_bit!(
        cloud_key::TFHECloudKey,
        result::TFHEEncryptedBit, a::TFHEEncryptedBit, b::TFHEEncryptedBit,
        lsb_carry::TFHEEncryptedBit, tmp::TFHEEncryptedBit)
    tfhe_gate_XNOR!(cloud_key, tmp, a, b)
    tfhe_gate_MUX!(cloud_key, result, tmp, lsb_carry, a)
end

# this function compares two multibit words, and puts the max in result
function encrypted_minimum!(
        cloud_key::TFHECloudKey, result::Array{TFHEEncryptedBit, 1},
        a::Array{TFHEEncryptedBit, 1}, b::Array{TFHEEncryptedBit, 1})

    nb_bits = length(result)

    params = tfhe_parameters(cloud_key)
    tmps = [TFHEEncryptedBit(params) for i in 1:2]

    # initialize the carry to 0
    tfhe_gate_CONSTANT!(cloud_key, tmps[1], Int32(0))
    # run the elementary comparator gate n times
    for i in 1:nb_bits
        encrypted_compare_bit!(cloud_key, tmps[1], a[i], b[i], tmps[1], tmps[2])
    end

    # tmps[1] is the result of the comparaison: 0 if a is larger, 1 if b is larger
    # select the max and copy it to the result
    for i in 1:nb_bits
        tfhe_gate_MUX!(cloud_key, result[i], tmps[1], b[i], a[i])
    end
end


function process(cloud_key, ciphertext1, ciphertext2)

    # if necessary, the params are inside the key
    params = tfhe_parameters(cloud_key)

    # do some operations on the ciphertexts: here, we will compute the
    # minimum of the two
    result = [TFHEEncryptedBit(params) for i in 1:16]
    encrypted_minimum!(cloud_key, result, ciphertext1, ciphertext2)

    result
end
=#


function process(cloud_key, ciphertext1, ciphertext2)
    params = tfhe_parameters(cloud_key)
    result = empty_ciphertext(params, size(ciphertext1)...)
    tfhe_gate_AND!(cloud_key, result, ciphertext1, ciphertext2)
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
