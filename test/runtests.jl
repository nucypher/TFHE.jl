using Base.Iterators: product

using Jute
using Random
using TFHE


gate_tests = [
    ("NAND", gate_nand, 2, !&),
    ("OR", gate_or, 2, |),
    ("AND", gate_and, 2, &),
    ("XOR", gate_xor, 2, xor),
    ("XNOR", gate_xnor, 2, (x, y) -> xor(x, ~y)),
    ("NOT", gate_not, 1, ~),
    ("NOR", gate_nor, 2, !|),
    ("ANDNY", gate_andny, 2, (x, y) -> (~x) & y),
    ("ANDYN", gate_andyn, 2, (x, y) -> x & (~y)),
    ("ORNY", gate_orny, 2, (x, y) -> (~x) | y),
    ("ORYN", gate_oryn, 2, (x, y) -> x | (~y)),
    ("MUX", gate_mux, 3, (x, y, z) -> x ? y : z),
]

gate_test_ids = [gate_test[1] for gate_test in gate_tests]


@testcase "gate" for gate_test in (gate_tests => gate_test_ids)
    rng = MersenneTwister(123)
    secret_key, cloud_key = make_key_pair(rng)

    _, gate, nargs, reference = gate_test

    for bits in product([(false, true) for i in 1:nargs]...)
        ebits = [encrypt(rng, secret_key, b) for b in bits]
        eres = gate(cloud_key, ebits...)
        res = decrypt(secret_key, eres)
        ref_res = reference(bits...)
        @test res == ref_res
    end

end


runtests()
