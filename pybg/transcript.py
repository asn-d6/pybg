"""
A basic Fiat-Shamir transcript object. It absorbs objects and spits out challenges.
"""

from py_ecc import optimized_bls12_381 as b
from bg_types import G1PointVector, FieldElementVector, FieldElement, G1Point
from hashlib import sha256

MODULUS = b.curve_order

def serialize_point(pt: G1Point):
    # Helper: Serializes an elliptic curve point.
    pt = b.normalize(pt)
    return pt[0].n.to_bytes(64, 'little') + pt[1].n.to_bytes(64, 'little')

def hash(x: bytes):
    return sha256(x).digest()

class Transcript:
    def __init__(self):
        self.digest = b""

    def absorb_points(self, ps: G1PointVector):
        """Add elliptic curve points to the transcript"""
        for p in ps:
            self.digest += serialize_point(p)

    def absorb_scalars(self, xs: FieldElementVector):
        """Add a bunch of scalars to the transcript"""
        for x in xs:
            self.digest += str(x).encode()

    def get_challenge_scalar(self) -> FieldElement:
        """Generate a scalar using the current state of the transcript"""
        challenge = int.from_bytes(hash(self.digest), 'little') % MODULUS
        # Add challenge to the digest. We do this so that we don't return the same challenge when this func is called
        # multiple times in a row
        self.digest += str(challenge).encode()
        return challenge
