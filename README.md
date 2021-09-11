# Prostate Shimming Project

## Background
Design shim coils to correct the off-resonance created by rectal air around the prostate region using MATLAB to provide a better resolution of prostate, which can improve prostate cancer identification. Some of the challenges are accounting for the anatomical difference between the physical phatom and real patients and the overcoming the extensive runtime of genetic algorithm (GA). Taking more patients' MRI scans into consideration and developing more parallel programming in GA code can lead to some potential improvement.

## Procedure
1. Load Field Maps (B0)
2. Define region of interest
3. Plot the prostate and the body
4. Generate Coil Coordinates
5. Calculate Bz with Biot Savart Law
6. Derive shim current
7. Shimming
8. Genetic algorithm

## Installation
1. Download MATLAB from [Go to the Support Web Site](https://support.west-wind.com)

## Results:
* Improve prostate field homogeneity by 60%

## function order:
1. Apply_Manual_Mask
2. Plot_Coils_9_vars
3. Biot_Savart
4. Solve_Shim_Current
5. Bz_Calc
