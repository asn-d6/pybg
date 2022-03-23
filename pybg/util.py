from py_ecc import optimized_bls12_381 as b
import random

MODULUS = b.curve_order

def msm(pts: list, scalars: list):
    """Naive linear combination of a list of points and values (aka multiscalar multiplication)"""
    assert len(pts) == len(scalars)

    o = b.Z1
    for pt, value in zip(pts, scalars):
        o = b.add(o, b.multiply(pt, value))
    return o

def is_power_of_two(x):
    return x and (x & (x-1) == 0)

def inv(a):
    """Modular inverse using eGCD algorithm"""
    if a == 0:
        return 0
    lm, hm = 1, 0
    low, high = a % MODULUS, MODULUS
    while low > 1:
        r = high // low
        nm, new = hm - lm * r, high - low * r
        lm, low, hm, high = nm, new, lm, low
    return lm % MODULUS

def get_inner_product(a, b):
    assert len(a) == len(b)
    return sum(x * y % MODULUS for x, y in zip(a, b)) % MODULUS

# Returns the (left|right) half of a container
def left_half(x):
    return x[:len(x)//2]
def right_half(x):
    return x[len(x)//2:]

def apply_permutation(a, perm):
    """Return permuted container `a` using the permutation `perm`"""
    assert len(a) == len(perm)
    return [a[i] for i in perm]
