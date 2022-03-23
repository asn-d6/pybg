"""
Grand-product argument prover
"""

import random

from py_ecc import optimized_bls12_381 as b

from bg_types import G1Point, FieldElement, G1PointVector, FieldElementVector
import inner_product_prove as ipa_prove
import gprod
from transcript import Transcript
from util import msm, get_inner_product, inv

MODULUS = b.curve_order

def prove(transcript: Transcript, crs_vec_G: G1PointVector, crs_U: G1Point,
          A: G1Point, gprod_result: FieldElement,
          vec_a: FieldElementVector, n_blinders: int) -> gprod.GrandProductProof:
    """
    Prove that there exists `vec_a` such that:
    - `A` is a commitment to `vec_a`
    - `gprod_result` is the product of the non-blinder elements of `vec_a`
    """
    n = len(crs_vec_G)
    ell = n - n_blinders

    # Step 1
    vec_b = []
    vec_b.append(1)
    for i, a in enumerate(vec_a[1:ell]):
        vec_b.append((a * vec_b[-1]) % MODULUS)
    vec_b.extend([random.randint(0, MODULUS) for _ in range(n_blinders)])

    B = msm(crs_vec_G, vec_b)
    bl = get_inner_product(vec_a[ell:], vec_b[ell:])

    transcript.absorb_points([A, B])
    transcript.absorb_scalars([bl])
    x = transcript.get_challenge_scalar()
    inv_x = inv(x)

    # Step 2
    # Start building C
    C = msm(crs_vec_G[:ell], [1]*ell)
    C = b.multiply(C, MODULUS - inv_x)
    C = b.add(C, A)

    vec_c = []
    pow_x = x
    pow_x2 = 1 # the second x in each c
    for a in vec_a[1:ell]: # iterate over a vector but skip the first element
        c = (a * pow_x - pow_x2) % MODULUS
        vec_c.append(c)
        pow_x = (pow_x * x) % MODULUS
        pow_x2 = (pow_x2 * x) % MODULUS

    # Append c with a_1 in the end
    c = (vec_a[0] * pow_x - pow_x2) % MODULUS
    vec_c.append(c)

    # Finally add blinders to the c vector
    pow_x = (pow_x * x) % MODULUS
    for a in vec_a[ell:]:
        vec_c.append((a * pow_x) % MODULUS)

    # Build the new basis
    crs_H = []
    pow_inv_x = inv_x
    for G in crs_vec_G[1:ell]:
        H = b.multiply(G, pow_inv_x)
        crs_H.append(H)
        pow_inv_x = (pow_inv_x * inv_x) % MODULUS
    crs_H.append(b.multiply(crs_vec_G[0], pow_inv_x))

    # Also add blinders to crs_H
    pow_inv_x = (pow_inv_x * inv_x) % MODULUS
    for G in crs_vec_G[ell:]:
        crs_H.append(b.multiply(G, pow_inv_x))

    # Step 3
    inner_prod = (bl * (x ** (ell+1)) + gprod_result * (x ** ell) - 1) % MODULUS
    ipa_proof = ipa_prove.prove(transcript, crs_vec_G, crs_H, crs_U, B, C, inner_prod, vec_b, vec_c)

    # Sanity check
    assert get_inner_product(vec_b, vec_c) == inner_prod

    return gprod.GrandProductProof(B, bl, ipa_proof)
