push!(LOAD_PATH, "../src")

using TFHE
using Random


function int_to_bits(x::T) where T <: Integer
    [((x >> (i - 1)) & 1 != 0) for i in 1:(sizeof(T)*8)]
end


function bits_to_int(::Type{T}, x::Array{Bool, 1}) where T <: Integer
    result = zero(T)
    for i in 1:min(sizeof(T) * 8, length(x))
        result |= (x[i]<<(i-1))
    end
    result
end


function encrypt()

    rng = MersenneTwister(123)
    secret_key, cloud_key = tfhe_key_pair(rng)

    # generate encrypt the 16 bits of 2017
    plaintext1 = UInt16(2017)
    bits1 = int_to_bits(plaintext1)
    ciphertext1 = [tfhe_encrypt_bit(rng, secret_key, bits1[i]) for i in 1:16]

    # generate encrypt the 16 bits of 42
    plaintext2 = UInt16(42)
    bits2 = int_to_bits(plaintext2)
    ciphertext2 = [tfhe_encrypt_bit(rng, secret_key, bits2[i]) for i in 1:16]

    secret_key, cloud_key, ciphertext1, ciphertext2
end


# elementary full comparator gate that is used to compare the i-th bit:
#   input: ai and bi the i-th bit of a and b
#          lsb_carry: the result of the comparison on the lowest bits
#   algo: if (a==b) return lsb_carry else return b
function encrypted_compare_bit(
        cloud_key::TFHECloudKey,
        a::TFHEEncryptedBit, b::TFHEEncryptedBit, lsb_carry::TFHEEncryptedBit)
    tmp = tfhe_gate_XNOR(cloud_key, a, b)
    tfhe_gate_MUX(cloud_key, tmp, lsb_carry, a)
end

# this function compares two multibit words, and puts the max in result
function encrypted_minimum(
        cloud_key::TFHECloudKey,
        a::Array{TFHEEncryptedBit, 1}, b::Array{TFHEEncryptedBit, 1})

    nb_bits = length(a)

    # initialize the carry to 0
    tmps1 = tfhe_gate_CONSTANT(cloud_key, false)
    # run the elementary comparator gate n times
    for i in 1:nb_bits
        tmps1 = encrypted_compare_bit(cloud_key, a[i], b[i], tmps1)
    end

    # tmps1 is the result of the comparaison: 0 if a is larger, 1 if b is larger
    # select the max and copy it to the result
    [tfhe_gate_MUX(cloud_key, tmps1, b[i], a[i]) for i in 1:nb_bits]
end


function process(cloud_key, ciphertext1, ciphertext2)
    # do some operations on the ciphertexts: here, we will compute the
    # minimum of the two
    encrypted_minimum(cloud_key, ciphertext1, ciphertext2)
end


function verify(secret_key, answer)
    # decrypt and rebuild the 16-bit plaintext answer
    bits = [tfhe_decrypt_bit(secret_key, answer[i]) for i in 1:length(answer)]
    int_answer = bits_to_int(UInt16, bits)

    println("Answer: $int_answer")
end


secret_key, cloud_key, ciphertext1, ciphertext2 = encrypt()
answer = process(cloud_key, ciphertext1, ciphertext2)
verify(secret_key, answer)
