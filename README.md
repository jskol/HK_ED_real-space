# HK_ED_real-space
A code to do an exact diagonalization calculations on a one-dimensional chain with Hatsugai-Kohmoto and Hubbard types of interaction
## Necessary libraries
This code uses following libraries 
1. libcommute for second quatization and Hilbertspace construction
     https://github.com/krivenko/libcommute
2. lambda-lanczos for Lanczos procedure of finding n-lowest eigenstates and eigenvalues
    https://github.com/mrcdr/lambda-lanczos

In addition, to speed-up the calculation support for the following is available
1. OpenMP (parallel programming) 