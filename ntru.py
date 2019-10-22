#!/usr/bin/env python3
"""NTRU v0.1

Usage:
  ntru.py [options] enc PUB_KEY_FILE [FILE]
  ntru.py [options] dec PRIV_KEY_FILE [FILE]
  ntru.py [options] gen N P Q PRIV_KEY_FILE PUB_KEY_FILE
  ntru.py [options] att PUB_KEY_FILE [FILE]
  ntru.py (-h | --help)
  ntru.py --version

Options:
  -h, --help         Show this screen.
  --version          Show version.
  -d, --debug        Debug mode.
  -v, --verbose      Verbose mode.

"""
from docopt import docopt
from mathutils import *
from sympy.abc import x
from sympy import ZZ, Poly
from padding import *
import numpy as np
import logging
import sys
import math
import os
from sympy.polys.polyerrors import NotInvertible
from collections import namedtuple

debug = False
verbose = False

class NTRU(NamedTuple):
    N: int
    p: int
    q: int

    def R(self):
        return Poly(x**self.N - 1, x).st_domain(ZZ)

def decrypt_method(cipher, f, f_p, q, p, R):
    """ Decrypts the message with the private key f """
    a_poly = ((f * cipher) % R).trunc(q)
    return ((f_p * a_poly) % R).trunc(p)


def generate_random_keys(N, p, q):
    """ Randomly generates the public and privates keys. Returns Numpy arrays"""
    R = Poly(x ** N - 1, x).set_domain(ZZ)
    g = random_poly(N, int(math.sqrt(N)))
    for i in range(10):
        f = random_poly(N, N // 3, neg_ones_diff=-1)
        try:
            f_p = invert_poly(f, R, p)
            f_q = invert_poly(f, R, q)
            break
        except NotInvertible as ex:
            pass
    h = ((p * f_q * g) % R).trunc(q)
    return asArr(h,N), asArr(f,N), asArr(f_p,N), asArr(g,N)


def generate_keys(N, p, q, priv_key_file, pub_key_file):
    """ Generates the public and private keys """
    # Generate the keys and convert them to numpy array form
    h, f, f_p, g = generate_random_keys(N, p, q)
    # Store the results in a file
    np.savez_compressed(priv_key_file, N=N, p=p, q=q, f=f, f_p=f_p, g=g)
    np.savez_compressed(pub_key_file, N=N, p=p, q=q, h=h)


def encrypt(pub_key_file, msg):
    """ Encrypt the message using the public key """
    # Load the NTRU system
    pub_key = np.load(pub_key_file, allow_pickle=True)
    N,p,q = int(pub_key['N']), int(pub_key['p']), int(pub_key['q'])
    R = Poly(x ** N - 1, x).set_domain(ZZ)
    h = asPoly(pub_key['h'].astype(np.int))

    # Pad the message and break up into blocks
    msg = padding_encode(msg, N).reshape((-1, N))
    def encrypt_method(msg_poly):
        """ Encrypts the message with the public key """
        blinding_poly = random_poly(N, int(math.sqrt(q)))
        return (((blinding_poly * h) + msg_poly) % R).trunc(q)
    output = np.concatenate([asArr(encrypt_method(asPoly(b)),N) for b in msg])

    output = [[0 if c == '0' else 1 for c in np.binary_repr(n, width=int(math.log2(q)))] for n in output]
    return np.array(output).flatten()


def asPoly(arr):
    return Poly(arr[::-1],x).set_domain(ZZ)


def asArr(poly, N):
    tmp = poly.all_coeffs()[::-1]
    return np.pad(tmp, (0, N - len(tmp)))


def decrypt(priv_key_file, cipher):
    """ Decrypt the message """
    priv_key = np.load(priv_key_file, allow_pickle=True)
    N,p,q = int(priv_key['N']), int(priv_key['p']), int(priv_key['q'])
    f = asPoly(priv_key['f'].astype(np.int))
    f_p = asPoly(priv_key['f_p'].astype(np.int))
    R = Poly(x ** N - 1, x).set_domain(ZZ)

    #return please_work(cipher, f, f_p, q, p, R, N)

    k = int(math.log2(q))
    pad = 0 if (len(cipher) % k) == 0 else k - (len(cipher) % k)
    cipher = np.array([int("".join(n.astype(str)), 2) for n in
                          np.pad(np.array(cipher), (0, pad), 'constant').reshape((-1, k))])
    cipher = cipher.reshape((-1, N))
    output = np.array([])
    for b in cipher:
        next_output = decrypt_method(asPoly(b), f, f_p, q, p, R).all_coeffs()[::-1]
        if len(next_output) < N:
            next_output = np.pad(next_output, (0, N - len(next_output)), 'constant')
        output = np.concatenate((output, next_output))
    return padding_decode(output, N)


def parseIntoProg(N, q, Lat):
    """ Write a magma program to solve the lattice """
    magProg = (f'print "{40 * "="}\\n";\n'
               f'print "NTRU Lattice\\n";\n;'
               f'print "{40 * "="}\\n";\n'
               f'B := MatrixRing( IntegerRing(), {2*N}) ! {Lat};\n'
               f'print "Unreduced NTRU Lattice";\n'
               f'print "NTRU Lattice after LLL";\n'
               f'Lattice(B);\n')

    with open("magma_Code.txt", 'w') as fileobj:
        fileobj.write(magProg)

def genNTRULattice(N, h, p, q):
    """ Generates the NTRU Lattice """

    inverse = lambda i: i if (p * i) % q == 1 else inverse(i+1)

    p_q = inverse(2) # Start searching for the inverse between 2 .. p
    # Top left quadrant has a diagonal of ones
    arr = np.eye(2*N, dtype=int)
    # Puts the value of q along the bottom right diagonal.
    for i in range(N, 2*N):
        arr[i, i] = q
    # The rotations of the public key in the top right quadrant
    for i in range(N):
        arr[i][N:] = np.hstack((h[-i:], h[:-i]))

    arr[:N, N:] *= p_q  # Multiply the top left quadrant by p_q
    return list(np.ndarray.flatten(arr))

def parseOutput(output):
    """ Parses the output from Magma and returns a Numpy array with each row
    being one of the reduced basis elements"""
    beginning = output.find('(')
    vectors_str = output.replace('\n',' ').replace('\t', ' ')
    vectors_arr = []
    while vectors_str.find('(') != -1:
        start = vectors_str.find('(')
        end = vectors_str.find(')')
        vectors_arr.append(vectors_str[start+1:end])
        vectors_str = vectors_str[end+1:]
    return np.array([np.array(v.split()[1:-1]).astype(int) for v in
                            vectors_arr])

def searchForPrivKey(latBasis, h, N, p, q, cipher):
    """ Searches for the private key in the reduced basis vectors """
    R = Poly(x ** N - 1, x).set_domain(ZZ)
    inverse = lambda i: i if (p * i) % q == 1 else inverse(i+1)
    p_q = inverse(2)
    for i, row in enumerate(latBasis, start=1):
        if abs(np.amax(row)) <= 1 and abs(np.sum(row)) <= 1:
            maybe_f = asPoly(row[:N])
            attempt = asArr((p_q * (maybe_f * h) % R).trunc(q), N)
            if abs(np.amax(attempt)) <= 1 and abs(np.sum(attempt)) <= 1:
                try:
                    Maybe_h, f_p = generate_public_key(maybe_f, asPoly(row[N:]), R, p, q)
                except NotInvertible:
                    print("Not Invertible")
                    continue
                #if Maybe_h == h:
                print("Attempt to decrypt")
                f_Arr = asArr(maybe_f, N)
                print(f"\nAttempting to decrypt with potential private key {f_Arr}\n")
                return please_work(cipher, maybe_f, f_p, q, p, R, N)
                #sys.stdout.buffer.write(np.packbits(np.array(yeet).astype(np.int)).tobytes())


def please_work(cipher, f, f_p, q, p, R, N):
    k = int(math.log2(q))
    pad = 0 if (len(cipher) % k) == 0 else k - (len(cipher) % k)
    cipher = np.array([int("".join(n.astype(str)), 2) for n in
                          np.pad(np.array(cipher), (0, pad), 'constant').reshape((-1, k))])
    cipher = cipher.reshape((-1, N))
    output = np.array([])
    for b in cipher:
        next_output = decrypt_method(asPoly(b), f, f_p, q, p, R).all_coeffs()[::-1]
        if len(next_output) < N:
            next_output = np.pad(next_output, (0, N - len(next_output)), 'constant')
        output = np.concatenate((output, next_output))
    return padding_decode(output, N)

def latticeAttack(pubKeyFile, cipher):
    """ Attacks the NTRU system via lattice basis reduction """
    # Read in the public data
    pub_key = np.load(pubKeyFile, allow_pickle=True)
    N, p, q = int(pub_key['N']), int(pub_key['p']), int(pub_key['q'])
    h = pub_key['h'].astype(np.int)

    # Generate the lattice
    ntruLat = genNTRULattice(N, h, p, q)
    parseIntoProg(N, q, ntruLat)

    # Call magma
    result = os.popen("magma -b < magma_Code.txt")
    bases = parseOutput(result.read())
    result.close()

    decrypted_msg = searchForPrivKey(bases, asPoly(h), N, p, q, cipher)


    priv_key = np.load('key_priv.npz', allow_pickle=True)
    f = priv_key['f'].astype(np.int)
    f_p = priv_key['f_p'].astype(np.int)
    g = priv_key['g'].astype(np.int)
    R = Poly(x ** N - 1, x).set_domain(ZZ)
    print("\nProper Decryption\n")
    msg = please_work(cipher, asPoly(f), asPoly(f_p), q, p, R, N)
    sys.stdout.buffer.write(np.packbits(np.array(msg).astype(np.int)).tobytes())
    #print(f"\nf = {list(f)}")
    #print(f"g = {list(g)}")
    #print(f"h = {h}")
    print("\nMy Decryption\n")
    return decrypted_msg


if __name__ == '__main__':
    #args = docopt(__doc__, version='NTRU v0.1')
    #root = logging.getLogger()
    #root.setLevel(logging.DEBUG)
    #ch = logging.StreamHandler(sys.stdout)
    #if args['--debug']:
    #    ch.setLevel(logging.DEBUG)
    #elif args['--verbose']:
    #    ch.setLevel(logging.INFO)
    #else:
    #    ch.setLevel(logging.WARN)
    #root.addHandler(ch)

    #log.debug(args)
    input_arr, output = None, None
    if not args['gen']:
        if args['FILE'] is None or args['FILE'] == '-':
            input = sys.stdin.read() if poly_input else sys.stdin.buffer.read()
        else:
            with open(args['FILE'], 'rb') as file:
                input = file.read()
        input_arr = np.unpackbits(np.frombuffer(input, dtype=np.uint8))
        input_arr = np.trim_zeros(input_arr, 'b')

    if args['gen']:

        generate_keys(int(args['N']), int(args['P']), int(args['Q']), args['PRIV_KEY_FILE'], args['PUB_KEY_FILE'])
    elif args['enc']:
        output = encrypt(args['PUB_KEY_FILE'], input_arr)
    elif args['dec']:
        output = decrypt(args['PRIV_KEY_FILE'], input_arr)
    elif args['att']:
        output = latticeAttack(args['PUB_KEY_FILE'], input_arr)

    if not args['gen']:
        sys.stdout.buffer.write(np.packbits(np.array(output).astype(np.int)).tobytes())

