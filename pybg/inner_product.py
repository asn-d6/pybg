"""
Inner product argument verifier
"""

from py_ecc import optimized_bls12_381 as b
from dataclasses import dataclass

from bg_types import G1Point, FieldElement, G1PointVector
from transcript import Transcript
from util import msm, inv, left_half, right_half

MODULUS = b.curve_order

@dataclass
class IPAProof():
    R: G1Point
    S: G1Point
    bl_1: FieldElement
    bl_2: FieldElement
    vec_B_L: G1PointVector
    vec_B_R: G1PointVector
    vec_C_L: G1PointVector
    vec_C_R: G1PointVector
    tip_b: FieldElement
    tip_c: FieldElement

def verify(transcript: Transcript, crs_vec_G: G1PointVector, crs_vec_H: G1PointVector, crs_U: G1Point,
           B: G1Point, C: G1Point, z: FieldElement, proof: IPAProof) -> bool:
    """
    Verify that `z` is the inner product of the vectors commited in `B` and `C`.
    """
    # Step 1
    transcript.absorb_points([B, C, proof.R, proof.S])
    transcript.absorb_scalars([z, proof.bl_1, proof.bl_2])
    x = transcript.get_challenge_scalar()

    z = (z + x*proof.bl_1 + (x**2)*proof.bl_2) % MODULUS
    B = b.add(B, b.multiply(proof.R, x))
    C = b.add(C, b.multiply(proof.S, x))

    # Step 2
    transcript.absorb_scalars([x])
    x = transcript.get_challenge_scalar()

    U = b.multiply(crs_U, x)
    B = b.add(B, b.multiply(U, z))

    # Step 3
    for i in range(len(proof.vec_B_L)):
        G_L, G_R = left_half(crs_vec_G), right_half(crs_vec_G)
        H_L, H_R = left_half(crs_vec_H), right_half(crs_vec_H)

        transcript.absorb_points([proof.vec_B_L[i], proof.vec_C_L[i], proof.vec_B_R[i], proof.vec_C_R[i]])
        x = transcript.get_challenge_scalar()
        x_inv = inv(x)

        B = msm([proof.vec_B_L[i], B, proof.vec_B_R[i]], [x, 1, x_inv])
        C = msm([proof.vec_C_L[i], C, proof.vec_C_R[i]], [x, 1, x_inv])

        crs_vec_G = [b.add(GL, b.multiply(GR, x_inv)) for (GL, GR) in zip(G_L, G_R)]
        crs_vec_H = [b.add(HL, b.multiply(HR, x)) for (HL, HR) in zip(H_L, H_R)]

    # Step 4
    assert len(crs_vec_G) == len(crs_vec_H) == 1
    exp_B = msm([crs_vec_G[0], U], [proof.tip_b, proof.tip_b * proof.tip_c])
    exp_C = b.multiply(crs_vec_H[0], proof.tip_c)

    assert b.eq(B, exp_B) and b.eq(C, exp_C)

    return True
