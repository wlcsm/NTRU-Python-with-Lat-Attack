# NTRU Implementation in Python 3

This is a basic implementation of the NTRU Cryptosystem completed as university coursework.

*Disclaimer*: This is not a complete or secure implementation and should not be
used as anything other than an educational tool.

Some of the details of this system were taken from https://github.com/jkrauze/ntru , such as certain polynomial coercion operations, and parts of the commandline interface. However, the majority of the program is completely different.

## Useage

```
./ntru.py gen N p q [PRIVATE_KEY_FILE] [PUBLIC_KEY_FILE]
./ntru.py enc [PUBLIC_KEY_FILE [CIPHERTEXT] [MESSAGE_FILE]
./ntru.py dec [PRIVATE_KEY_FILE] [CIPHERTEXT_FILE]
./ntru.py att [PUBLIC_KEY_FILE] [CIPHERTEXT_FILE]
```
The best explanation of how to use the system is to look at the testing script outlined below.

## Sample usage

The bash script "test.sh" generates and encrypts a message contained in the "message.txt" file.
It then decrypts the ciphertext using the public key, then attacks the cryptosystem with a lattice attack.

Though it should be able to handle any size of text, the program is quite slow and a large plaintext could take a while to encrypt and decrypt.

The attack is only particularly successful for N less that 45.


William Cashman 2019
