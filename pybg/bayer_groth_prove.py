"""
Prover of the Bayer-Groth shuffle argument
"""

import math, random

from py_ecc import optimized_bls12_381 as b

from bayer_groth import ShuffleCRS, ShuffleProof
from bg_types import FieldElement, G1PointVector
import gprod_prove, sameexp_prove, multiexp_prove
from util import msm, get_inner_product, apply_permutation
from transcript import Transcript

MODULUS = b.curve_order

# Number of blinders we need in the shuffle proof
N_BLINDERS = 4

def prove(crs: ShuffleCRS,
          vec_R: G1PointVector, vec_S: G1PointVector, vec_T: G1PointVector, vec_U: G1PointVector,
          permutation: list, r: FieldElement) -> ShuffleProof:
    """
    Proves that there exist `permutation` and `r` such that:

    The elements of `vec_R` and `vec_S` were permuted using `permutation` and randomized using `r`, and the results are
    in `vec_T` and `vec_U` respectively.
    """
    # Number of non-blinder elements used in this proof
    ell = len(vec_R)
    # Total number of elements used in proof (including blinders)
    n = N_BLINDERS + ell

    transcript = Transcript() # Our Fiat-Shamir transcript

    # Step 1
    vec_s_blinders = [random.randint(0, MODULUS) for _ in range(N_BLINDERS)]
    vec_perm_with_s_blinders = permutation + vec_s_blinders
    M = msm(crs.vec_G, vec_perm_with_s_blinders)

    transcript.absorb_points(vec_T + vec_U + [M])
    vec_a = [transcript.get_challenge_scalar() for _ in range(ell)]

    # Step 2
    # Add a bunch of blinders to `a` vector
    vec_a_blinders = [random.randint(0, MODULUS) for _ in range(N_BLINDERS)]
    vec_a_permuted_with_blinders = apply_permutation(vec_a, permutation) + vec_a_blinders

    A = msm(crs.vec_G, vec_a_permuted_with_blinders)

    transcript.absorb_points([A])
    alpha, beta = transcript.get_challenge_scalar(), transcript.get_challenge_scalar()

    # Step 3
    # We use `vec_perm_with_s_blinders` here so that the blinders follow the permuted numbers
    permuted_polynomial_factors = [(a + m * alpha + beta) % MODULUS for a,m in zip(vec_a_permuted_with_blinders, vec_perm_with_s_blinders)]
    # We compute the grand product over the non-blinder part of the polynomial factors
    gprod_result = math.prod(permuted_polynomial_factors[:ell]) % MODULUS
    A_1 = msm([A, M] + crs.vec_G, [1, alpha] + [beta]*n)
    gprod_proof = gprod_prove.prove(transcript, crs.vec_G, crs.U, A_1, gprod_result, permuted_polynomial_factors, N_BLINDERS)

    # Sanity check: make sure that permuted polynomial has same roots as the regular polynomial
    # vec_a_with_blinders = vec_a + vec_a_blinders
    # polynomial_factors = [a + m*alpha + beta for a,m in zip(vec_a_with_blinders, list(range(ELL)) + vec_s_blinders)]
    # assert gprod_result == (math.prod(permuted_polynomial_factors[:ELL]) % MODULUS)

    # Step 4
    transcript.absorb_points([A])
    vec_gamma, vec_delta = [], [] # need...more...blinders
    for _ in range(N_BLINDERS):
        vec_gamma.append(transcript.get_challenge_scalar())
        vec_delta.append(transcript.get_challenge_scalar())

    R = msm(vec_R, vec_a)
    S = msm(vec_S, vec_a)
    r_t = get_inner_product(vec_gamma, vec_a_blinders)
    r_u = get_inner_product(vec_delta, vec_a_blinders)
    T = msm([R, crs.G_t], [r, r_t])
    U = msm([S, crs.G_u], [r, r_u])

    sameexp_proof = sameexp_prove.prove(transcript, crs.G_t, crs.G_u, R, S, T, U, r, r_t, r_u)

    # Step 5
    vec_T_with_blinders = vec_T + [b.multiply(crs.G_t, gamma) for gamma in vec_gamma]
    vec_U_with_blinders = vec_U + [b.multiply(crs.G_u, delta) for delta in vec_delta]
    multiexp_proof = multiexp_prove.prove(transcript, crs.vec_G, vec_T_with_blinders, vec_U_with_blinders, A, T, U, vec_a_permuted_with_blinders)

    return ShuffleProof(M, A, T, U, gprod_proof, sameexp_proof, multiexp_proof)
