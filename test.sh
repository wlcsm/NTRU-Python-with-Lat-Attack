#!/bin/bash
# N = 23 works 25 maybe works
./ntru.py gen 103 3 64 key_priv key_pub
echo "Generated keys"
./ntru.py enc key_pub.npz message.txt > ciphertext.txt
echo "Encrypted message"
./ntru.py dec key_priv.npz ciphertext.txt
echo "Decrypted message"

