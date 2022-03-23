"""
Same-exponentiation argument verifier
"""

from py_ecc import optimized_bls12_381 as b
from dataclasses import dataclass

from bg_types import G1Point, FieldElement
from transcript import Transcript
from util import msm

MODULUS = b.curve_order

@dataclass
class SameExponentProof():
    B_t: G1Point
    B_u: G1Point
    z_r: FieldElement
    z_t: FieldElement
    z_u: FieldElement

def verify(transcript: Transcript, crs_G_t: G1Point, crs_G_u: G1Point,
           R: G1Point, S: G1Point, T: G1Point, U: G1Point, proof: SameExponentProof) -> bool:
    """
    Verify proof that there exists `r`, `r_t` and `r_u` s.t.:
    - `T = r * R + r_t * G_t`
    - `U = r * S + r_u * G_u`
    """
    # Step 1
    transcript.absorb_points([R, S, T, U])
    transcript.absorb_scalars([proof.B_t, proof.B_u])
    x = transcript.get_challenge_scalar()

    # Step 2
    expected_1 = msm([proof.B_t, T, R, crs_G_t], [1, x, MODULUS - proof.z_r, MODULUS - proof.z_t])
    expected_2 = msm([proof.B_u, U, S, crs_G_u], [1, x, MODULUS - proof.z_r, MODULUS - proof.z_u])
    assert b.eq(expected_1, b.Z1) and b.eq(expected_2, b.Z1)

    return True
