#!/bin/bash
# The main factor which cause decryption error is the value of q.
# Through testing I have found q = 64 to be a stable number for N < 64.
# If one starts seeming decryption errors simply increase the value of q.
echo "Generating Keys"
./ntru.py gen 43 3 64 key_priv key_pub
echo "Encrypting message"
./ntru.py enc key_pub.npz ciphertext message.txt
echo "Decrypting message"
./ntru.py dec key_priv.npz ciphertext.npz
echo "Attacking NTRU"
./ntru.py att key_pub.npz ciphertext.npz

