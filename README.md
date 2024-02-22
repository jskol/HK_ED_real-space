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

## How-to quick guide
The code in the current form allows to calculate multiple quntites using proper flags. These quatities and their flags in paranthesis:
1. Single particle spectral function (-A)
2. Local spin-excitation spectrum (-S)
3. Two point correlator (-E)
4. Spin-spin correlator (-E)
In addition, following modifications to the Hamiltonian can be switched on
1. change from H-K to Hubbard interaction (-H)
2. Change from Open Boundary conditions to Periodic Boundary Conditions (-P)

The implemented Hamiltonian parameters are :
1. interaction strength U (-U)
2. hopping t (-t)
3. second hopping tp, needed for SSH-type models (-p). By defult it will be the same as t
4. chainlenght N (-N)
5. Linear Voltage across the chain (-V)
6. Zeeman Magnetic field (-M) 

Example:
To set up calulations for a N=8 chain with HK interaction U=5 hopping t=1 and calulating site-resolved spectral function use
./HK-Edge_states -U 5 -N 8 -A -t -1 -p -1
The results will be in "finite_system_N_8_t_-1_tp_-1_U_5.dat"
