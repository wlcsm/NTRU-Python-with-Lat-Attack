#!/usr/bin/env python3
"""NTRU v0.1

Usage:
  ntru.py gen N P Q PRIV_KEY_FILE PUB_KEY_FILE
  ntru.py enc PUB_KEY_FILE CIPHERTEXT_FILE FILE
  ntru.py dec PRIV_KEY_FILE FILE
  ntru.py att PUB_KEY_FILE FILE
"""

from docopt import docopt
from mathutils import *
from ntruPolyOps import *
from sympy import ZZ, Poly
from sympy.abc import x
from sympy.polys.polyerrors import NotInvertible
import numpy as np
import sys
import math
import os


def generate_keys(ntru, priv_key_file, pub_key_file):
    """ Generates the public and private keys """
    # Generate the keys and convert them to numpy array form
    N = ntru.N
    p = ntru.p
    q = ntru.q

    g = random_poly(N, int(math.sqrt(N)))
    for i in range(10):
        f = random_poly(N, N // 3, neg_ones_diff=-1)
        try:
            f_p = invert_poly(f, ntru.R, p)
            f_q = invert_poly(f, ntru.R, q)
            break
        except NotInvertible as ex:
            pass
    h = ((p * f_q * g) % ntru.R).trunc(q)

    h, f, f_p, g = asArr(h,N), asArr(f,N), asArr(f_p,N), asArr(g,N)

    # Store the results in a file
    np.savez_compressed(priv_key_file, N=N, p=p, q=q, f=f, f_p=f_p, g=g)
    np.savez_compressed(pub_key_file, N=N, p=p, q=q, h=h)


def encrypt(pub_key_file, msg):
    """ Encrypt the message using the public key """
    ntru, h = load_NTRU(pub_key_file, isPriv=False)

    # Pad the message and break up into blocks
    padding = np.zeros(2*ntru.N - (len(msg) % ntru.N))
    padding[ntru.N + 1] = 1
    msg_blocks = np.concatenate((msg, padding)).reshape((-1, ntru.N))

    # Encrypts the message with the public key
    def enc_poly(msg_poly):
        e = random_poly(ntru.N, int(math.sqrt(ntru.q)))
        return (((e * h) + msg_poly) % ntru.R).trunc(ntru.q)

    # Encrypt all the blocks and recombine into one array
    return np.concatenate([asArr(enc_poly(asPoly(b)),ntru.N) for b in msg_blocks])


def decrypt(priv_key_file, cipher):
    """ Decrypts the ciphertext """
    ntru, f, f_p = load_NTRU(priv_key_file, isPriv=True)

    # Organise into blocks
    cipher_blocks = cipher.reshape((-1, ntru.N))

    def dec_poly(cipher):
        """ Decrypts the message with the private key """
        f_m = ((f * cipher) % ntru.R).trunc(ntru.q)
        return ((f_p * f_m) % ntru.R).trunc(ntru.p)

    output = np.concatenate([asArr(dec_poly(asPoly(b)), ntru.N) for b in cipher_blocks])

    find_last_zero = lambda n: n if output[n] == 0 else find_last_zero(n - 1)
    return output[:find_last_zero(len(output)-1)]


def parseIntoProg(N, q, Lat):
    """ Write a magma program to solve the lattice """
    magProg = (f'print "{40 * "="}\\n";\n'
               f'print "NTRU Lattice\\n";\n;'
               f'print "{40 * "="}\\n";\n'
               f'B := MatrixRing( IntegerRing(), {2*N}) ! {Lat};\n'
               f'print "Unreduced NTRU Lattice";\n'
               f'print "NTRU Lattice after LLL";\n'
               f'V := Lattice(B);\n'
               f'LLL(V:Proof:=false, Delta:=0.9999999, Eta:=0.5000001);')

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

    arr[:N, N:] *= p_q  # Multiply the top right quadrant by p_q
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
    print(max(latBasis[0]))
    R = Poly(x ** N - 1, x).set_domain(ZZ)
    inverse = lambda i: i if (p * i) % q == 1 else inverse(i+1)
    p_q = inverse(2)
    for i, row in enumerate(latBasis, start=1):
        if abs(np.amax(row)) <= 1 and abs(np.sum(row)) <= 1:
            maybe_f = asPoly(row[:N])
            attempt = asArr((p_q * (maybe_f * h) % R).trunc(q), N)
            if abs(np.amax(attempt)) <= 1 and abs(np.sum(attempt)) <= 1:
                try:
                    f_q = invert_poly(maybe_f, R, q)
                    f_p = invert_poly(maybe_f, R, p)
                    h = (p * f_q * asPoly(row[N:]) % R).trunc(q)
                except NotInvertible:
                    print("Not Invertible")
                    continue
                #if Maybe_h == h:
                print("Attempt to decrypt")
                f_Arr = asArr(maybe_f, N)
                print(f"\nAttempting to decrypt with potential private key {f_Arr}\n")
                return please_work(cipher, maybe_f, f_p, q, p, R, N)
                #sys.stdout.buffer.write(np.packbits(np.array(yeet).astype(np.int)).tobytes())
    return "Not found"


def please_work(cipher, f, f_p, q, p, R, N):
    cipher = cipher.reshape((-1, N))

    def decrypt_method(cipher):
        """ Decrypts the message with the private key f """
        a_poly = ((f * cipher) % R).trunc(q)
        return ((f_p * a_poly) % R).trunc(p)


    output = np.concatenate(([asArr(decrypt_method(asPoly(b)), N) for b in cipher]))

    find_last_zero = lambda n: n if output[n] == 0 else find_last_zero(n - 1)
    return output[:find_last_zero(len(output)-1)]


def latticeAttack(pubKeyFile, cipher):
    """ Attacks the NTRU system via lattice basis reduction """
    # Read in the public data
    ntru, h = load_NTRU(pubKeyFile, isPriv=False)

    # Generate the lattice
    ntruLat = genNTRULattice(ntru.N, asArr(h,ntru.N), ntru.p, ntru.q)
    parseIntoProg(ntru.N, ntru.q, ntruLat)

    # Call magma
    result = os.popen("magma -b < magma_Code.txt")
    bases = parseOutput(result.read())
    result.close()
    print(bases[0])
    print([sum([x*x for x in row]) for row in bases])

    decrypted_msg = searchForPrivKey(bases, h, ntru.N, ntru.p, ntru.q, cipher)
    if decrypted_mg == "Not Found":
        return "Not found"

    print(f"Yeet {decrypted_msg}")
    priv_key = np.load('key_priv.npz', allow_pickle=True)
    f = priv_key['f'].astype(np.int)
    print(sum([x*x for x in f]))
    f_p = priv_key['f_p'].astype(np.int)
    g = priv_key['g'].astype(np.int)
    print(sum([x*x for x in g]))
    print("\nProper Decryption\n")
    msg = please_work(cipher, asPoly(f), asPoly(f_p), ntru.q, ntru.p, ntru.R, ntru.N)
    sys.stdout.buffer.write(np.packbits(np.array(msg).astype(np.int)).tobytes())
    #print(f"\nf = {list(f)}")
    #print(f"g = {list(g)}")
    #print(f"h = {h}")
    print("\nMy Decryption\n")
    return decrypted_msg


if __name__ == '__main__':
    args = docopt(__doc__, version='NTRU v0.1')
    input_arr, output = None, None

    if args['enc']:
        if args['FILE'] is None or args['FILE'] == '-':
            input = sys.stdin.read() if poly_input else sys.stdin.buffer.read()
        else:
            with open(args['FILE'], 'rb') as file:
                input = file.read()
        input_arr = np.unpackbits(np.frombuffer(input, dtype=np.uint8))
        input_arr = np.trim_zeros(input_arr, 'b')

    if args['gen']:
        NTRU_Par = NTRU(int(args['N']), int(args['P']), int(args['Q']))
        generate_keys(NTRU_Par, args['PRIV_KEY_FILE'], args['PUB_KEY_FILE'])
    elif args['enc']:
        output = encrypt(args['PUB_KEY_FILE'], input_arr)
        np.savez_compressed(args['CIPHERTEXT_FILE'], cipher=output)
    elif args['dec']:
        input_arr = np.load(args['FILE'], allow_pickle=True)['cipher'].astype(int)
        output = decrypt(args['PRIV_KEY_FILE'], input_arr)
    elif args['att']:
        input_arr = np.load(args['FILE'], allow_pickle=True)['cipher'].astype(int)
        output = latticeAttack(args['PUB_KEY_FILE'], input_arr)

    if not args['gen'] and not args['enc']:
        sys.stdout.buffer.write(np.packbits(np.array(output).astype(np.int)).tobytes())

