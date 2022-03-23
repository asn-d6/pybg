"""
Grand-product argument verifier
"""

from py_ecc import optimized_bls12_381 as b
from dataclasses import dataclass

import inner_product as ipa
from bg_types import G1Point, FieldElement, G1PointVector
from transcript import Transcript
from util import msm, inv

MODULUS = b.curve_order

@dataclass
class GrandProductProof():
    B: G1Point
    bl: FieldElement
    ipa_proof: ipa.IPAProof

def verify(transcript: Transcript, crs_vec_G: G1PointVector, crs_U: G1Point,
           A: G1Point, gprod_result: FieldElement, n_blinders: int, proof: GrandProductProof) -> bool:
    """
    Verify that `gprod_result` is the product of the non-blinder vector elements commited in `A`.
    """
    n = len(crs_vec_G)
    ell = n - n_blinders

    # Step 1
    transcript.absorb_points([A, proof.B])
    transcript.absorb_scalars([proof.bl])
    x = transcript.get_challenge_scalar()
    inv_x = inv(x)

    # Step 2
    # Start building C
    C = msm(crs_vec_G[:ell], [1]*ell)
    C = b.multiply(C, MODULUS - inv_x)
    C = b.add(C, A)

    # Now build the new basis
    crs_H = []
    pow_inv_x = inv_x
    for G in crs_vec_G[1:ell]:
        H = b.multiply(G, pow_inv_x)
        crs_H.append(H)
        pow_inv_x = pow_inv_x * inv_x % MODULUS
    crs_H.append(b.multiply(crs_vec_G[0], pow_inv_x))

    # Also add blinders to crs_H
    pow_inv_x = pow_inv_x * inv_x % MODULUS
    for G in crs_vec_G[ell:]:
        crs_H.append(b.multiply(G, pow_inv_x))

    # Step 3
    inner_prod = (proof.bl * (x ** (ell+1)) + gprod_result * (x ** ell) - 1) % MODULUS
    return ipa.verify(transcript, crs_vec_G, crs_H, crs_U, proof.B, C, inner_prod, proof.ipa_proof)

