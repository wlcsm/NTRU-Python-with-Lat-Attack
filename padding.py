""" Padding functions for the NTRU system.

    Padding is as follows: Given a message and a block size, the system will pad
    the message with zeros until it can be broken up into blocks. Then it will
    attach another block on the end. If the original message was padded with k
    zeros, then the last block will contain k ones at the end of the block, and
    the rest are zeros.
"""
import numpy as np

def padding_encode(msg, block_size):
    """ Break up the input arr into blocks of size 'block_size' with padding
    with zeros. It then adds a block onto the end to indicate what the padding
    was."""
    # The stuff we need to pad on the end
    padding = ([0] * block_size) + ([1] * (block_size - len(msg) % block_size))
    return np.concatenate((msg, np.array(padding)))


def padding_decode(cipher, block_size):
    """ Removes the padding on the cipher """
    find_last_zero = lambda n: n if cipher[n] == 0 else find_last_zero(n - 1)
    return cipher[:find_last_zero(-1)]
