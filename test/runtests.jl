using Base.Iterators: product

using Jute
using Random
using TFHE


gate_tests = [
    ("NAND", tfhe_gate_NAND, 2, !&),
    ("OR", tfhe_gate_OR, 2, |),
    ("AND", tfhe_gate_AND, 2, &),
    ("XOR", tfhe_gate_XOR, 2, xor),
    ("XNOR", tfhe_gate_XNOR, 2, (x, y) -> xor(x, ~y)),
    ("NOT", tfhe_gate_NOT, 1, ~),
    ("NOR", tfhe_gate_NOR, 2, !|),
    ("ANDNY", tfhe_gate_ANDNY, 2, (x, y) -> (~x) & y),
    ("ANDYN", tfhe_gate_ANDYN, 2, (x, y) -> x & (~y)),
    ("ORNY", tfhe_gate_ORNY, 2, (x, y) -> (~x) | y),
    ("ORYN", tfhe_gate_ORYN, 2, (x, y) -> x | (~y)),
    ("MUX", tfhe_gate_MUX, 3, (x, y, z) -> x ? y : z),
]

gate_test_ids = [gate_test[1] for gate_test in gate_tests]


@testcase "gate" for gate_test in (gate_tests => gate_test_ids)
    rng = MersenneTwister(123)
    secret_key, cloud_key = tfhe_key_pair(rng)

    _, gate, nargs, reference = gate_test

    for bits in product([(false, true) for i in 1:nargs]...)
        ebits = [tfhe_encrypt_bit(rng, secret_key, b) for b in bits]
        eres = gate(cloud_key, ebits...)
        res = tfhe_decrypt_bit(secret_key, eres)
        ref_res = reference(bits...)
        @test res == ref_res
    end

end


runtests()
