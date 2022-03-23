"""
Verifier for a Bayer-Groth shuffle argument
"""

import math

from py_ecc import optimized_bls12_381 as b
from dataclasses import dataclass

from bg_types import G1Point, G1PointVector
import gprod, sameexp, multiexp
from util import msm
from transcript import Transcript

MODULUS = b.curve_order

# Number of blinders we need in the shuffle proof
N_BLINDERS = 4

@dataclass
class ShuffleCRS:
    """
    The CRS of the shuffle proof. Includes basis elements used in the various subarguments.
    """
    vec_G: G1PointVector
    U: G1Point
    G_t: G1Point
    G_u: G1Point


@dataclass
class ShuffleProof:
    M: G1Point
    A: G1Point
    T: G1Point
    U: G1Point
    gprod_proof: gprod.GrandProductProof
    sameexp_proof: sameexp.SameExponentProof
    multiexp_proof: multiexp.MultiExpProof


def verify(crs: ShuffleCRS,
           vec_R: G1PointVector, vec_S: G1PointVector, vec_T: G1PointVector, vec_U: G1PointVector, proof: ShuffleProof) -> bool:
    """
    Verifies that the elements of `vec_R` and `vec_S` were permuted and randomized, and
    the output is in `vec_T` and `vec_U` respectively.
    """
    # Number of non-blinder elements used in this proof
    ell = len(vec_R)
    # Total number of elements used in proof (including blinders)
    n = ell + N_BLINDERS

    # Get our Fiat-Shamir transcript
    transcript = Transcript()

    # Step 1
    transcript.absorb_points(vec_T + vec_U + [proof.M])
    vec_a = [transcript.get_challenge_scalar() for _ in range(ell)]

    # Step 2
    transcript.absorb_points([proof.A])
    alpha, beta = transcript.get_challenge_scalar(), transcript.get_challenge_scalar()

    # Step 3
    polynomial_coeffs = [(a + i * alpha + beta) % MODULUS for i,a in enumerate(vec_a)]
    gprod_result = math.prod(polynomial_coeffs) % MODULUS
    A_1 = msm([proof.A, proof.M] + crs.vec_G, [1, alpha] + [beta]*n)
    assert gprod.verify(transcript, crs.vec_G, crs.U, A_1, gprod_result, N_BLINDERS, proof.gprod_proof)

    # Step 4
    transcript.absorb_points([proof.A])
    vec_gamma, vec_delta = [], [] # need...more...blinders
    for _ in range(N_BLINDERS):
        vec_gamma.append(transcript.get_challenge_scalar())
        vec_delta.append(transcript.get_challenge_scalar())

    R = msm(vec_R, vec_a)
    S = msm(vec_S, vec_a)
    assert sameexp.verify(transcript, crs.G_t, crs.G_u, R, S, proof.T, proof.U, proof.sameexp_proof)

    # Step 5
    vec_T_with_blinders = vec_T + [b.multiply(crs.G_t, gamma) for gamma in vec_gamma]
    vec_U_with_blinders = vec_U + [b.multiply(crs.G_u, delta) for delta in vec_delta]
    assert multiexp.verify(transcript, crs.vec_G, vec_T_with_blinders, vec_U_with_blinders, proof.A, proof.T, proof.U, proof.multiexp_proof)

    return True

