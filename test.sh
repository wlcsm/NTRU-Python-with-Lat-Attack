#!/bin/bash
# N = 23 works 25 maybe works
./ntru.py gen 43 3 64 key_priv key_pub
echo "Generated keys"
./ntru.py enc key_pub.npz ciphertext message.txt
echo "Encrypted message"
./ntru.py dec key_priv.npz ciphertext.npz
echo "Decrypted message"

