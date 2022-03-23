"""
Multi-exponentiation argument verifier
"""

from py_ecc import optimized_bls12_381 as b
from dataclasses import dataclass

from bg_types import G1Point, FieldElement, G1PointVector
from transcript import Transcript
from util import msm, inv, left_half, right_half

MODULUS = b.curve_order

@dataclass
class MultiExpProof():
    R: G1Point
    T_bl: G1Point
    U_bl: G1Point
    vec_T_L: G1PointVector
    vec_T_R: G1PointVector
    vec_U_L: G1PointVector
    vec_U_R: G1PointVector
    vec_C_L: G1PointVector
    vec_C_R: G1PointVector
    tip_a: FieldElement

def verify(transcript: Transcript, crs_G: G1PointVector,
           vec_T: G1PointVector, vec_U: G1PointVector, A: G1Point, T: G1Point, U: G1Point, proof: MultiExpProof) -> bool:
    """
    Verify that:
    - `A` is a commitment to vector `vec_a`
    - `T` is the result of an MSM between `vec_T` and `vec_a`
    - `U` is the result of an MSM between `vec_U` and `vec_a`
    """
    # Step 1
    transcript.absorb_points([A, T, U, proof.R, proof.T_bl, proof.U_bl])
    x = transcript.get_challenge_scalar()

    A = b.add(A, b.multiply(proof.R, x))
    T = b.add(T, b.multiply(proof.T_bl, x))
    U = b.add(U, b.multiply(proof.U_bl, x))

    # Step 2: log(n) rounds of recursion
    for i in range(len(proof.vec_C_L)):
        T_L, T_R = left_half(vec_T), right_half(vec_T)
        U_L, U_R = left_half(vec_U), right_half(vec_U)
        G_L, G_R = left_half(crs_G), right_half(crs_G)

        transcript.absorb_points([proof.vec_T_L[i], proof.vec_U_L[i], proof.vec_T_R[i],
                                  proof.vec_U_R[i], proof.vec_C_L[i], proof.vec_C_R[i]])
        x = transcript.get_challenge_scalar()
        x_inv = inv(x)

        A = msm([proof.vec_C_L[i], A, proof.vec_C_R[i]], [x, 1, x_inv])
        T = msm([proof.vec_T_L[i], T, proof.vec_T_R[i]], [x, 1, x_inv])
        U = msm([proof.vec_U_L[i], U, proof.vec_U_R[i]], [x, 1, x_inv])

        crs_G = [b.add(GL, b.multiply(GR, x)) for (GL, GR) in zip(G_L, G_R)]
        vec_T = [b.add(TL, b.multiply(TR, x)) for (TL, TR) in zip(T_L, T_R)]
        vec_U = [b.add(UL, b.multiply(UR, x)) for (UL, UR) in zip(U_L, U_R)]

    # Step 3
    assert len(crs_G) == len(vec_T) == len(vec_U) == 1
    exp_A = b.multiply(crs_G[0], proof.tip_a)
    exp_T = b.multiply(vec_T[0], proof.tip_a)
    exp_U = b.multiply(vec_U[0], proof.tip_a)
    assert b.eq(A, exp_A) and b.eq(T, exp_T) and b.eq(U, exp_U)

    return True

