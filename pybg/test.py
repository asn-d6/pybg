"""
End-to-end tests for all the zero-knowledge arguments involved.
"""

import unittest, time
import random
import math

from py_ecc import optimized_bls12_381 as b

import gprod, sameexp, multiexp, bayer_groth, inner_product as ipa
import gprod_prove, sameexp_prove, multiexp_prove, bayer_groth_prove, inner_product_prove as ipa_prove
from transcript import Transcript
from util import get_inner_product, apply_permutation, msm

MODULUS = b.curve_order

def gen_generator_points(count):
    """
    Create `count` generator points for BLS12-381
    """
    BLS12_381_COFACTOR = 76329603384216526031706109802092473003
    points = []
    x = b.FQ(1)
    while len(points) < count:
        y = (x ** 3 + b.b) ** ((b.field_modulus + 1) // 4)
        if b.is_on_curve((x, y, b.FQ(1)), b.b):
            points.append(b.multiply((x, y, b.FQ(1)), BLS12_381_COFACTOR))
        x += b.FQ(1)
    return points

def get_random_permutation(n):
    """Return a random permutation with `n` elements"""
    p = list(range(n))
    random.shuffle(p)
    return p

time_cache = [time.time()]
def get_time_delta():
    time_cache.append(time.time())
    return time_cache[-1] - time_cache[-2]

# Number of total group elements involved in the shuffle proof (including zero-knowledge blinders)
N = 128
# Number of blinder elements needed for the shuffle proof
N_BLINDERS = 4
# Number of actual useful non-blinder elements involved in the shuffle proof
ELL = N - N_BLINDERS

class TestInnerProductArgument(unittest.TestCase):
    def test_inner_product_argument(self):
        generators = gen_generator_points(2*N + 1)
        crs_G = generators[:N]
        crs_H = generators[N:-1]
        crs_U = generators[-1]
        print("ipa: generated generator points: {:.3f}s".format(get_time_delta()))

        # Create dummy vectors to use in this test
        vec_b = [random.randint(0, MODULUS) for i in range(N)]
        vec_c = [random.randint(0, MODULUS) for i in range(N)]
        z = get_inner_product(vec_b, vec_c)

        B = msm(crs_G, vec_b)
        C = msm(crs_H, vec_c)
        print("ipa: generated random vectors and their commitments: {:.3f}s".format(get_time_delta()))

        # Create ZKP that vec_b \dot vec_c == z
        proof = ipa_prove.prove(Transcript(), crs_G[:N], crs_H[:N], crs_U, B, C, z, vec_b, vec_c)
        print("ipa: proof generated: {:.3f}s".format(get_time_delta()))

        assert ipa.verify(Transcript(), crs_G, crs_H, crs_U, B, C, z, proof)
        print("ipa: proof verified: {:.3f}s".format(get_time_delta()))

class TestMultiExpProof(unittest.TestCase):
    def test_multi_exp_argument(self):
        # Create generators needed for multiexp proof
        generators = gen_generator_points(3*N)
        crs_G = generators[:N]
        vec_T = generators[N:2*N]
        vec_U = generators[2*N:]
        print("multiexp: generated generator points: {:.3f}s".format(get_time_delta()))

        # Create witnesses
        vec_a = [random.randint(0, MODULUS) for i in range(N)]
        # Commitment to witnesses
        A = msm(crs_G, vec_a) # commitment to A

        # Create commitments (multiexp results) that verifier knows
        T = msm(vec_T, vec_a)
        U = msm(vec_U, vec_a)
        print("multiexp: generated witness vector and its commitments: {:.3f}s".format(get_time_delta()))

        get_time_delta()
        # Create multiexp ZKP
        proof = multiexp_prove.prove(Transcript(), crs_G, vec_T, vec_U, A, T, U, vec_a)
        print("multiexp: proof generated: {:.3f}s".format(get_time_delta()))

        assert multiexp.verify(Transcript(), crs_G, vec_T, vec_U, A, T, U, proof)
        print("multiexp: proof verified: {:.3f}s".format(get_time_delta()))

class TestGrandProduct(unittest.TestCase):
    def test_grand_product_argument(self):
        # Create generators
        generators = gen_generator_points(N+1) # make some extra generator points for the blinders
        crs_G = generators[:N]
        crs_U = generators[-1]
        print("gprod: generated generator points: {:.3f}s".format(get_time_delta()))

        # Create dummy `a` vector to use in this test
        vec_a = [random.randint(0, MODULUS) for i in range(ELL)]
        # Compute the grand product with the current vector (before adding the blinders)
        gprod_result = math.prod(vec_a) % MODULUS
        # Add a bunch of blinders to `vec_a`
        vec_a.extend([random.randint(0, MODULUS) for _ in range(N_BLINDERS)])

        # Compute commitment to `vec_a`
        A_1 = msm(crs_G, vec_a)
        print("gprod: generated random vector, blinders and its commitments: {:.3f}s".format(get_time_delta()))

        # Create grand product proof
        gprod_proof = gprod_prove.prove(Transcript(), crs_G, crs_U, A_1, gprod_result, vec_a, N_BLINDERS)
        print("gprod: proof generated: {:.3f}s".format(get_time_delta()))

        # Verify the proof
        assert gprod.verify(Transcript(), crs_G, crs_U, A_1, gprod_result, N_BLINDERS, gprod_proof)
        print("gprod: proof verified: {:.3f}s".format(get_time_delta()))

class SameExpProof(unittest.TestCase):
    def test_same_exp_argument(self):
        # Create generators
        generators = gen_generator_points(4)
        crs_G_t = generators[0]
        crs_G_u = generators[1]
        print("sameexp: generated generator points: {:.3f}s".format(get_time_delta()))

        r = random.randint(0, MODULUS)
        r_t = random.randint(0, MODULUS)
        r_u = random.randint(0, MODULUS)

        R = generators[2]
        T = msm([R, crs_G_t], [r, r_t])

        S = generators[3]
        U = msm([S, crs_G_u], [r, r_u])

        print("sameexp: generated scalars and group elements: {:.3f}s".format(get_time_delta()))

        sameexp_proof = sameexp_prove.prove(Transcript(), crs_G_t, crs_G_u, R, S, T, U, r, r_t, r_u)
        print("sameexp: proof generated: {:.3f}s".format(get_time_delta()))

        assert sameexp.verify(Transcript(), crs_G_t, crs_G_u, R, S, T, U, sameexp_proof)
        print("sameexp: proof verified: {:.3f}s".format(get_time_delta()))


class TestShuffleProof(unittest.TestCase):
    def test_shuffle_argument(self):
        """
        Run an end-to-end test of the Bayer-Groth argument.
        Generate some input; shuffle it; get the proof; verify it.
        """
        # Create the basis of our CRS
        generators = gen_generator_points(5*N)
        print("bg: generated generator points: {:.3f}s".format(get_time_delta()))

        crs_G = generators[:N]
        crs_U = generators[-1]
        crs_G_t = generators[-2]
        crs_G_u = generators[-3]
        crs = bayer_groth.ShuffleCRS(crs_G, crs_U, crs_G_t, crs_G_u)

        # Generate a random permutation and the randomizing factor
        permutation = get_random_permutation(ELL)
        r = random.randint(0, MODULUS) # randomizing factor

        # Create the input vectors
        vec_R = generators[N:N+ELL]
        vec_S = generators[N+ELL:N+2*ELL]
        assert(len(vec_R) == len(vec_S) == ELL)

        # Create the output vectors by randomizing and permuting the input vectors
        vec_T = [b.multiply(R_i, r) for R_i in vec_R]
        vec_T = apply_permutation(vec_T, permutation)
        vec_U = [b.multiply(S_i, r) for S_i in vec_S]
        vec_U = apply_permutation(vec_U, permutation)
        print("bg: finished shuffling and randomizing: {:.3f}s".format(get_time_delta()))

        # Create a shuffle proof and verify it
        shuffle_proof = bayer_groth_prove.prove(crs, vec_R, vec_S, vec_T, vec_U, permutation, r)
        print("bg: finished shuffle proof: {:.3f}s".format(get_time_delta()))

        assert bayer_groth.verify(crs, vec_R, vec_S, vec_T, vec_U, shuffle_proof)
        print("bg: finished verifying shuffle proof: {:.3f}s".format(get_time_delta()))

if __name__ == '__main__':
    unittest.main()

