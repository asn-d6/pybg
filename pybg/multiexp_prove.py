"""
Multi-exponentiation argument prover
"""

from py_ecc import optimized_bls12_381 as b

import random

import multiexp
from bg_types import G1Point, G1PointVector, FieldElementVector
from transcript import Transcript
from util import msm, is_power_of_two, inv, left_half, right_half

MODULUS = b.curve_order

def prove(transcript: Transcript, crs_G: G1PointVector,
          vec_T: G1PointVector, vec_U: G1PointVector, A: G1Point, T: G1Point, U: G1Point,
          vec_a: FieldElementVector) -> multiexp.MultiExpProof:
    """
    Prove that there exists `vec_a` such that:
    - `A` is a commitment to `vec_a`
    - `T` is the result of an MSM between `vec_T` and `vec_a`
    - `U` is the result of an MSM between `vec_U` and `vec_a`
    """
    n = len(crs_G)
    assert len(crs_G) == len(vec_T) == len(vec_U) == len(vec_a)
    assert is_power_of_two(len(vec_a))

    vec_T_L, vec_T_R, vec_U_L, vec_U_R, vec_C_L, vec_C_R = [], [], [], [], [], []

    # Step 1
    vec_r = [random.randint(0, MODULUS) for i in range(n)]
    R = msm(crs_G, vec_r)
    T_bl = msm(vec_T, vec_r)
    U_bl = msm(vec_U, vec_r)

    transcript.absorb_points([A, T, U, R, T_bl, U_bl])
    x = transcript.get_challenge_scalar()

    # Rewrite the vectors b and c
    vec_a = [(a + x * r) % MODULUS for (a, r) in zip(vec_a, vec_r)]

    # Step 2: log(n) rounds of recursion
    while len(vec_a) > 1:
        a_L, a_R = left_half(vec_a), right_half(vec_a)
        T_L, T_R = left_half(vec_T), right_half(vec_T)
        U_L, U_R = left_half(vec_U), right_half(vec_U)
        G_L, G_R = left_half(crs_G), right_half(crs_G)

        Z_L_T = msm(T_R, a_L)
        Z_L_U = msm(U_R, a_L)
        Z_R_T = msm(T_L, a_R)
        Z_R_U = msm(U_L, a_R)

        C_L = msm(G_R, a_L)
        C_R = msm(G_L, a_R)

        # Append to proof
        vec_T_L.append(Z_L_T)
        vec_T_R.append(Z_R_T)
        vec_U_L.append(Z_L_U)
        vec_U_R.append(Z_R_U)
        vec_C_L.append(C_L)
        vec_C_R.append(C_R)

        transcript.absorb_points([Z_L_T, Z_L_U, Z_R_T, Z_R_U, C_L, C_R])
        x = transcript.get_challenge_scalar()
        x_inv = inv(x)

        # Generate half-size polynomial and points for the next round
        vec_a = [(aL + aR * x_inv) % MODULUS for (aL, aR) in zip(a_L, a_R)]
        vec_T = [b.add(TL, b.multiply(TR, x)) for (TL, TR) in zip(T_L, T_R)]
        vec_U = [b.add(UL, b.multiply(UR, x)) for (UL, UR) in zip(U_L, U_R)]
        crs_G = [b.add(GL, b.multiply(GR, x)) for (GL, GR) in zip(G_L, G_R)]

    # Step 3
    assert len(vec_a) == 1
    return multiexp.MultiExpProof(R, T_bl, U_bl, vec_T_L, vec_T_R, vec_U_L, vec_U_R, vec_C_L, vec_C_R, vec_a[0])
