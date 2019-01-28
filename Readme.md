# Numerical Linear Algebra: Numerical solutions to Stokes Equation

This repository implements the final project for the course Numerical Linear Algebra in Autumn 2018. Submitted version.



## Description of codes

The main codes of my implementations are stored under 'Matlab_Code'. Certain auxiliary programs which help us build the matrix for the problem are catogorized into `components`. The notations are consistent with my report. 



The main program is `VCycle.m` . It implements V-Cycle multigrid method to solve the saddle point problem, which is the discretization of the Stokes Equation. 

The smoothers contain:

- DGS-p (in Long Chen's paper): `DGS.m`.
- DGS-s (in the notes): `DGS_seq.m` and `DGS_seqM.m`. Their inner loop programs are `DGS_InnerC2.c` and `DGS_InnerM.m` . The former is implemented in C and the latter is implemenented in Matlab. `DGS_InnerC.c` is the history version of `DGS_InnerC2.c`. 
- Uzawa: `Uzawa.m`. We simply use \ in Matlab to solve the subproblem.
- Inexact Uzawa: `Inexact_Uzawa.m`. We call `CG_solver.m` to solve the subproblem.



We also implement Inexact Uzawa based on V-Cycle multi-grid. The main program is `ieuzawaVC.m`. We use V-Cycle multi-grid method to solve the subproblem. This is implemented in `VCycle_Inner.m`. The smoother for V-Cycle is Gauss-Seidel Iteration, which is implemented in `GS.m`.



In Appendix, we propose a modified V-Cycle multi-grid method based on DGS. The main program is `VCycle_mod.m`. The smoothers contain:

- DGS-p(in Long Chen's paper): `DGS_mod.m`.

- DGS-s (in the notes): `DGS_seqM_mod.m`. Its innner loop program is `DGS_InnerM_mod.m`



## How to run the codes

To run the codes, first, please use MEX in Matlab to compile `DGS_innerC2.c` in the following way:

`mex DGS_innerC2.c`

To recover the results for DGS-s and DGS-p, please run:

- $N=64, 128, 256$: `Test_DGS1.m` (DGS-s), `Test_DGS2.m` (DGS-p);
- $N=512, 1024, 2048$: `Test_DGS3.m` (DGS-s), `Test_DGS4.m` (DGS-p).

To recover the results for Uzawa, please run:

- $N = 64, 128$ : `Test_uzawa.m`

- $N = 256, 512$ : `Test_uzawa2.m`

- $N = 1024, 2048$ : `Test_uzawa3.m`

To recover the results for Inexact Uzawa, please run:

- $N = 64, 128$ : `Test_ieuzawa.m`

- $N = 256, 512$ : `Test_ieuzawa2.m`

- $N = 1024, 2048$ : `Test_ieuzawa3.m`

To recover the results for Inexact Uzawa based on V-Cycle, please run:

- $N = 64, 128$ : `Test_ieuzawaVC.m`

- $N = 256, 512$ : `Test_ieuzawaVC2.m`

- $N = 1024, 2048$ : `Test_ieuzawaVC3.m`



## Catogory

Our Matlab programs are generally catogorized as follows:

### bnd

'comp_*.m' generate the vector to describe the Neurman boundary condition

### comp

'comp_*.m' generate matrix which is component in the formulation of saddle point problem.

### spp

'comp_*.m' generate matrix in the saddle point problem

### test 

'test_*.m' implements the numerical experiements in the report

### true

'true_*.m' return the true value of the function

