# Nep2rat

'Nep2rat' is a tool to rationally approximate nonlinear matrix-valued functions, mainly aimed to solve nonlinear eigenvalue problems.

# Contents

This repository contains the function "nep2rat" which approximates the nonlinear functions, and all the experiments proposed in the referenced paper. Those experiments **require** the installation of the [NLEVP library](https://github.com/ftisseur/nlevp).

# Basic usage

     [Am, Bm, Rm] = nep2rat(F, Z)

returns the rational approximation Rm of the function F on the target set Z, and the linearized pencil (Am, Bm). The set Z must be given as a vector of points, while F could either be a function_handle or a struct, with fields:
    
  - F.coeffs: the matrix coefficients of F;
  - F.fun: the scalar functions that define f;

This structure mirrors the output of the NLEVP library. The default behaviour changes with the input F: if F is a function_handle, then nep2rat uses the "surrogate AAA with cyclic Leja--Bagby refinement" algorithm, while if F is a struct, it uses the "weighted AAA" algorithm. If the split form of F is available, we suggest to use the struct form of F as the input.

    [Am, Bm, Rm, info] = nep2rat(F, Z, opts) 
allows the user to specify many optional parameters, such as the precision of the approximation and the algorithms utilised. The output info is returns additional information, while opts is a structure with all the. An in-depth guide is contained in the help of nep2rat.

# References
Güttel, S., Negri Porzio, G.M. and Tisseur, F., 2020. [Robust rational approximations of nonlinear eigenvalue problems](http://eprints.maths.manchester.ac.uk/2796/1/gnt20.pdf). MIMS Eprints 2020.24
