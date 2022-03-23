"""
Inner product argument prover
"""

from py_ecc import optimized_bls12_381 as b

import random

import inner_product as ipa
from bg_types import G1Point, FieldElement, G1PointVector, FieldElementVector
from transcript import Transcript
from util import msm, is_power_of_two, inv, get_inner_product, left_half, right_half

MODULUS = b.curve_order

def prove(transcript: Transcript, crs_vec_G: G1PointVector, crs_vec_H: G1PointVector, crs_U: G1Point,
          B: G1Point, C: G1Point, z: FieldElement,
          vec_b: FieldElementVector, vec_c: FieldElementVector) -> ipa.IPAProof:
    """
    Prove that there exist `vec_b` and `vec_c` such that:
    - z is the inner product of `vec_b` and `vec_c`
    - B is the commitment of `vec_b`
    - C is the commitment of `vec_c`
    """
    n = len(vec_b)
    assert len(vec_b) == len(vec_c) == len(crs_vec_G) == len(crs_vec_H)
    assert is_power_of_two(len(vec_b))

    vec_B_L, vec_B_R, vec_C_L, vec_C_R = [], [], [], []

    # Step 1
    vec_r = [random.randint(0, MODULUS) for i in range(n)]
    vec_s = [random.randint(0, MODULUS) for i in range(n)]
    R = msm(crs_vec_G, vec_r)
    S = msm(crs_vec_H, vec_s)

    # Create blinders
    bl_1 = get_inner_product(vec_b, vec_s) + get_inner_product(vec_c, vec_r)
    bl_2 = get_inner_product(vec_r, vec_s)

    transcript.absorb_points([B, C, R, S])
    transcript.absorb_scalars([z, bl_1, bl_2])
    x = transcript.get_challenge_scalar()

    # Rewrite the vectors b and c
    vec_b = [(b + x*r) % MODULUS for (b, r) in zip(vec_b, vec_r)]
    vec_c = [(c + x*s) % MODULUS for (c, s) in zip(vec_c, vec_s)]

    # Step 2
    transcript.absorb_scalars([x])
    x = transcript.get_challenge_scalar()
    U = b.multiply(crs_U, x)

    # Step 3: log(n) rounds of recursion
    while len(vec_b) > 1:
        # Generate the left-side and right-side points
        b_L, b_R = left_half(vec_b), right_half(vec_b)
        c_L, c_R = left_half(vec_c), right_half(vec_c)
        G_L, G_R = left_half(crs_vec_G), right_half(crs_vec_G)
        H_L, H_R = left_half(crs_vec_H), right_half(crs_vec_H)

        C_L_b = b.add(msm(G_L, b_R), b.multiply(U, get_inner_product(b_R, c_L)))
        C_R_b = b.add(msm(G_R, b_L), b.multiply(U, get_inner_product(b_L, c_R)))
        C_L_c = msm(H_R, c_L)
        C_R_c = msm(H_L, c_R)

        # Append to proof
        vec_B_L.append(C_L_b)
        vec_C_L.append(C_L_c)
        vec_B_R.append(C_R_b)
        vec_C_R.append(C_R_c)

        transcript.absorb_points([C_L_b, C_L_c, C_R_b, C_R_c])
        x = transcript.get_challenge_scalar()
        x_inv = inv(x)

        vec_b = [(bL + bR * x) % MODULUS for (bL, bR) in zip(b_L, b_R)]
        vec_c = [(cL + cR * x_inv) % MODULUS for (cL, cR) in zip(c_L, c_R)]
        crs_vec_G = [b.add(GL, b.multiply(GR, x_inv)) for (GL, GR) in zip(G_L, G_R)]
        crs_vec_H = [b.add(HL, b.multiply(HR, x)) for (HL, HR) in zip(H_L, H_R)]

    # Step 4
    assert len(vec_b) == len(vec_c) == 1
    return ipa.IPAProof(R, S, bl_1, bl_2, vec_B_L, vec_B_R, vec_C_L, vec_C_R, vec_b[0], vec_c[0])

