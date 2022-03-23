# pybg (WIP)

**This is work in progress**

pybg is an implementation of the [Bayer-Groth zero-knowledge argument](http://www0.cs.ucl.ac.uk/staff/J.Groth/MinimalShuffle.pdf) for the correctness of a shuffle.

More specifically, pybg implements a [variant of the Bayer-Groth argument](https://github.com/ethresearch/Shuffle_SSLE/blob/master/docs/shuffle_ssle.pdf) designed by Mary Maller, which uses techniques from the Bulletproofs IPA to make the protocol simpler and more efficient.

## Purpose

This library is optimized for readability: it uses the same notation and structure as the [paper](https://github.com/ethresearch/Shuffle_SSLE/blob/master/docs/shuffle_ssle.pdf).

Our plan is to use the Bayer-Groth shuffle argument for a [leader election protocol](https://ethresear.ch/t/whisk-a-practical-shuffle-based-ssle-protocol-for-ethereum/11763#proofs-of-correct-shuffle-13) in Ethereum. pybg serves as a reference implementation of the shuffle argument, as well as a way to extract test vectors for other more optimized implementations of the Bayer-Groth argument.

pybg is not production-level software and is not optimized for performance.

pybg runs the argument over the BLS12-381 elliptic curve using the py_ecc library, which is what the Ethereum specs also use.

## Usage/Examples

The core of the shuffle argument can be found in `pybg/bayer_groth.py` and `pybg/bayer_groth_prove.py`.

See `test_shuffle_proof()` in `pybg/test.py` for a tutorial on how to use pybg's Bayer-Groth argument.

## Installation

You will need the `py_ecc` library to run pybg:

```bash
pip install py_ecc
```

## Running tests

To run the tests, run the following command:

```bash
    python pybg/test.py`
```
