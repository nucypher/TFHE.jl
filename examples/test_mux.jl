push!(LOAD_PATH, "../src")

using TFHE

if Base.thisminor(VERSION) > v"0.6"
    using Random
end


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


function reference_mux(bits1, bits2, bits3)
    BitArray(map((a, b, c) -> a ? b : c, bits1, bits2, bits3))
end


function encrypt()

    rng = MersenneTwister(123)

    tfhe_key_pair(rng)
    println("Key generation:")
    @time secret_key, cloud_key = tfhe_key_pair(rng)

    plaintext1 = Int16(2017)
    bits1 = int_to_bitarray(plaintext1)

    plaintext2 = Int16(42)
    bits2 = int_to_bitarray(plaintext2)

    plaintext3 = Int16(12345)
    bits3 = int_to_bitarray(plaintext3)

    ciphertext1 = tfhe_encrypt(rng, secret_key, bits1)
    ciphertext2 = tfhe_encrypt(rng, secret_key, bits2)
    ciphertext3 = tfhe_encrypt(rng, secret_key, bits3)

    println("Encryption:")
    ciphertext1, ciphertext2, ciphertext3 = @time begin
        (tfhe_encrypt(rng, secret_key, bits1),
        tfhe_encrypt(rng, secret_key, bits2),
        tfhe_encrypt(rng, secret_key, bits3))
    end

    secret_key, cloud_key, ciphertext1, ciphertext2, ciphertext3
end


function process(cloud_key, ciphertext1, ciphertext2, ciphertext3)
    params = tfhe_parameters(cloud_key)
    result = empty_ciphertext(params, size(ciphertext1)...)

    tfhe_gate_MUX!(cloud_key, result, ciphertext1, ciphertext2, ciphertext3)

    println("Processing:")
    @time tfhe_gate_MUX!(cloud_key, result, ciphertext1, ciphertext2, ciphertext3)

    result
end


function verify(secret_key, answer)
    tfhe_decrypt(secret_key, answer)
    println("Decryption:")
    @time answer_bits = tfhe_decrypt(secret_key, answer)

    int_answer = bitarray_to_int(answer_bits)
    println("Answer: $int_answer")
end


secret_key, cloud_key, ciphertext1, ciphertext2, ciphertext3 = encrypt()
answer = process(cloud_key, ciphertext1, ciphertext2, ciphertext3)
verify(secret_key, answer)
