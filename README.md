# HPS-SS-model

Coarse-Grained Dihedral angle potential for IDPs

CG simulations were performed using the LAMMPS molecular dynamics simulations package (Oct 2020 version), in which HPS-SS codes have been implemented.

Cα-based helix assignment rules
1) Helix fraction script (HelixFracDihed15.f90) is used to test different helix assignment rules
2) The output is then compared with the DSSP-based helix fraction of Ala40 peptide from atomistic simulation (AA-A40 folder)

Parameterization
1) A40 simulations with different **eps_d** values to figure out the reference value for highest helicity to match the experimental helical propensity for Alanine.
2) Helix propensity code was used to compute helix propensity (w) directly from the simulation trajectory.
3) **eps_d** for remaining 19 residues determined using a host-guest system such as A20XA20 or A20X4A20

Validation
1) TDP-43 CTD files was added as an example for the validation step
2) 
