push!(LOAD_PATH, "../src")

using TFHE
using Jute


key_pairs = @global_fixture begin

    rng = MersenneTwister(124)
    secret_key, cloud_key = tfhe_key_pair(rng)

    @produce (secret_key, cloud_key) "test_key"
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

    int_answer
end



@testcase "Key creation performance" begin

    rng = MersenneTwister(123)

    tfhe_key_pair(rng)

    tic()
    secret_key, cloud_key = tfhe_key_pair(rng)
    t = toq()

    # On iMac
    reference_time = 0.46
    @test_result (t / reference_time)
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


@testcase "MUX gate performance" for key_pair in key_pairs

    secret_key, cloud_key = key_pair
    params = secret_key.params

    bits1 = int_to_bitarray(Int16(2017))
    bits2 = int_to_bitarray(Int16(42))
    bits3 = int_to_bitarray(Int16(12345))
    rng = MersenneTwister(125)

    # Encryption

    tfhe_encrypt(rng, secret_key, bits1)
    tic()
    ciphertext1 = tfhe_encrypt(rng, secret_key, bits1)
    ciphertext2 = tfhe_encrypt(rng, secret_key, bits2)
    ciphertext3 = tfhe_encrypt(rng, secret_key, bits3)
    t = toq()

    # On iMac
    reference_time_enc = 0.00035
    @test_result (t / reference_time_enc)


    # Processing

    params = tfhe_parameters(cloud_key)
    result = empty_ciphertext(params, size(ciphertext1)...)

    tfhe_gate_MUX!(cloud_key, result, ciphertext1, ciphertext2, ciphertext3)

    tic()
    tfhe_gate_MUX!(cloud_key, result, ciphertext1, ciphertext2, ciphertext3)
    t = toq()

    # On iMac
    reference_time_process = 0.94
    @test_result (t / reference_time_process)


    # Decryption

    tfhe_decrypt(secret_key, result)
    tic()
    answer_bits = tfhe_decrypt(secret_key, result)
    t = toq()

    # On iMac
    reference_time_dec = 9e-6
    @test_result (t / reference_time_dec)

    int_answer = bitarray_to_int(answer_bits)
    ref_answer = bitarray_to_int(reference_mux(bits1, bits2, bits3))

    @test int_answer == ref_answer

end


runtests()
