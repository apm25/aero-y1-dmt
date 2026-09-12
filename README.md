# Wind Turbine BEM Solver (Imperial Aero Y1 DMT)

This repository contains a MATLAB script that I wrote, in collaboration with my team members for the Imperial College Aeronautics Year 1 DMT project. We had to design a wind turbine for a 10m/s wind tunnel with a 50cm max diameter. We achieved a peak power of 12W, with zero blade structural failure.

I built this solver to test different Design Tip Speed Ratios (TSRs) and automatically generate the spanwise chord and twist distributions. 

Note: The group's final presentation has been uploaded (`DMT_CTM_G19.pdf`) to show actual CAD renders, the manufacturing process, and how it performed in the wind tunnel.

## How the code works
It runs Blade Element Momentum theory across a range of wind speeds and TSRs to find the peak Cp. 
* It blends three cross-sections: a thick cylinder at the root (for structural integrity), an SG6040 airfoil for the mid-span, and an S834 airfoil at the tip.
* It dynamically interpolates the airfoil polar data based on the local Reynolds number.
* **Aero corrections:** It includes Prandtl tip/root losses, the Glauert empirical correction for high induction states, 3D stall delay near the hub, and Viterna extrapolation so the solver doesn't crash post-stall.

## How to use it
1. Keep `BEM.m` and the two `.csv` polar data files in the same folder.
2. Open `BEM.m` and adjust your rotor parameters at the top if you need to (it's currently hardcoded to R = 0.243m to safely clear the 50cm rig limit).
3. Run the script

It will print the physical chord constraints to the console and output a few `CAD_Blueprint_Blade_TSR*.csv` files for the given TSRs. These coordinates can be imported directly into SolidWorks or Fusion 360 to loft the 3D blade.
