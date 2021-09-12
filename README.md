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
1. Download and install [MATLAB](https://www.mathworks.com/login?uri=%2Fdownloads%2Fweb_downloads)
2. Install Image Processing Toolbox and Parallel Computing Toolbox under **Home** tab, **Add-ons**, **Get Add-ons**.
3. Save the five MATLAB functions to a folder and open the folder with MATLAB

## Usage
1. Supply your own phase maps to serve as B0 and magnitude maps for drawing mask
2. Use Apply_Manual_Mask to draw a mask around the region of interest (ROI). 
3. Come up with coil parameter, each coil should have 9 variables (specified in the header of the function)
    * the last variable, current, should be set to 1 for now
5. Feed the coil matrix into Plot_Coils_9_vars, observe coil position with respect to the ROI and body
6. Use BiotSavart with coil position, mask, and B0
7. Use Solve_Shim_Current with the output from BiotSavart and mask
8. Use Bz_Calc with the output from Solve_Shim_Current to get Bz

The functions are applied in the following order 
1. Apply_Manual_Mask
2. Plot_Coils_9_vars
3. Biot_Savart
4. Solve_Shim_Current
5. Bz_Calc

## Credits
Thanks to Hsin-Jung (Randy) Yang and Yu-Heng (Chris) Huang for coding structure and guidance 

## Testing input
The magnitude map **mag** and phase map **phasemap** can be found in the [testing physical phantom input](https://github.com/erictang0220/Cedars-Sinai-Internship/blob/main/UNIC_B0MapInVIVO_BH.mat)
