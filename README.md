# HPS-SS-model

Coarse-Grained Dihedral angle potential for IDPs

CG simulations were performed using the LAMMPS molecular dynamics simulations package (Oct 2020 version), in which HPS-SS codes have been implemented.

Cα-based helix assignment rules
1) Helix fraction script (HelixFracDihed15.f90) is used to test different helix assignment rules
The script requires xtc library to compile:
ifort -o gethelixfrac HelixFracDihed15.f90 -lxdrf -L/xtc

To run the script:
for (i,i+2) rule: ./gethelixfrac -x traj.xtc -o helix.dat -dihed 0.25 1.3 -assign 0110 -nn 1 -block 5 -eq 1000
for (i,i+4) rule: ./gethelixfrac -x traj.xtc -o helix.dat -dihed 0.25 1.7 -assign 0110 -nn 3 010 -block 5 -eq 1000

3) The output is then compared with the DSSP-based helix fraction of Ala40 peptide from atomistic simulation (AA-A40 folder)

Parameterization
1) A40 simulations with different **eps_d** values to figure out the reference value for highest helicity to match the experimental helical propensity for Alanine.
2) Helix propensity code (HelixPropensity.f90) was used to compute helix propensity (w) directly from the simulation trajectory.
to compile: ifort -o gethelixprop HelixPropensity.f90 -lxdrf -L/xtc
to run: ./gethelixprop -x alanine.xtc -o helicity.dat -e 10000 -t 0.25 1.7 -Nblock 5

4) **eps_d** for remaining 19 residues determined using a host-guest system such as A20XA20 or A20X4A20

Validation
1) TDP-43 CTD files was added as an example for the validation step
2) 
