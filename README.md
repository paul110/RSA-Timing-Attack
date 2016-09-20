# RSA-Timing-Attack

RSA attack based on execution time of the target machine

## Running the program

python attack.py ./target ./configuration

## DESCRIPTION

The file "configuration" represents a text file which contains the RSA modulus on the first line and the RSA public exponent on the second one, both represented as hexadecimal integer strings.

The file "target" is a binary which simulates a machine which takes a ciphertext as input and will decrypt it using a fixed unknown RSA key and return the corresponding plaintext as a result. At the same time the machine returns the time ( in clock cycles ) took for decryption.  

The program takes this as a parameter and  partially follows the paper  http://www.cs.jhu.edu/~fabian/courses/CS600.624/Timing-full.pdf in order to decipher the key used to decrypt he ciphertexts. It also implements some safety mechanisms in order to make the whole process more efficient, such as :

  * statistical threshold for guessing key bits

  * on bits with come up with warning the program simulates a few steps ahead in order to check whether it looks like it made the right choice on a certain bit of the key or not

  * program sets a stable key ( bits which are correct with very high probability ) in order to have a back up to return to if the key discovery looks like it is going to fail or does indeed fail

  * every time the process is backtracked, the program sends another bunch of ciphertexts to be decrypted in order to improve the statistical results when comparing the times. This also allows for the possibility that the program will recover the key using very few interactions, but accounts for the case in which it is not able to and it needs more.
