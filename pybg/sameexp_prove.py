"""
Same-exponentiation argument prover
"""

import random

from py_ecc import optimized_bls12_381 as b

import sameexp
from bg_types import G1Point, FieldElement
from transcript import Transcript
from util import msm

MODULUS = b.curve_order

def prove(transcript: Transcript, crs_G_t: G1Point, crs_G_u: G1Point,
          R: G1Point, S: G1Point, T: G1Point, U: G1Point,
          r: FieldElement, r_t: FieldElement, r_u: FieldElement) -> sameexp.SameExponentProof:
    """
    Prove that there exist `r`, `r_t` and `r_u` such that:
    - `T = r * R + r_t * G_t`
    - `U = r * S + r_u * G_u`
    """
    # Step 1
    bl_r = random.randint(0, MODULUS)
    bl_t = random.randint(0, MODULUS)
    bl_u = random.randint(0, MODULUS)

    B_t = msm([R, crs_G_t], [bl_r, bl_t])
    B_u = msm([S, crs_G_u], [bl_r, bl_u])

    transcript.absorb_points([R, S, T, U])
    transcript.absorb_scalars([B_t, B_u])
    x = transcript.get_challenge_scalar()

    # Step 2
    z_r = (bl_r + r * x) % MODULUS
    z_t = (bl_t + r_t * x) % MODULUS
    z_u = (bl_u + r_u * x) % MODULUS

    return sameexp.SameExponentProof(B_t, B_u, z_r, z_t, z_u)
