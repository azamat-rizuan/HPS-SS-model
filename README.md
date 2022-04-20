# HPS-SS-model

## Developing transferable bonded potentials for a coarse-grained model of Intrinsically Disordered Proteins

**Authors:** Azamat Rizuan [1], Nina Jovic [1], Young C. Kim [2], Jeetain Mittal [1]  
[1] Artie McFerrin Department of Chemical Engineering, Texas A&M University, College Station, Texas  
[2] Center for Materials Physics and Technology, Naval Research Laboratory, Washington, District of Columbia  

## Source codes:

CG simulations were performed using the LAMMPS molecular dynamics simulations package (Oct 2020 version), in which HPS-SS codes have been implemented. To run these codes successfully one would require LAMMPS Oct 2020 package installed with the files within 1.LAMMPS_subroutines.


### Input parameters of the HPS-SS model:
* The short-range van der Waals (vdW) and long-range electrostatic interactions between residues i and j is the same as HPS-Urry model (for more details https://doi.org/10.1002/pro.4094)
* Aminoacid hydropathy scale (aminoacids_hydropathy.dat)
* Van der Waals diameter of aminoacids (aminoacids_vdwdiameter.dat)
* Harmonic bond potential parameters
```
bond_style  harmonic
bond_coeff          1   10.000000    3.800000
```
* Angle potential parameters
```
angle_style bch
angle_coeff         1    4.300000
```
* Dihedral angle parameters (**eps_d**) for (i,i+2)  and (i,i+4) rules (eps_d_i+1_i+2.txt, eps_d_i_i+4.txt)
```
dihedral_style gaussian
dihedral_coeff        1    eps_d
```
**eps_d** for a given dihedral angle φ(i,i+3) is determined by the mixing rule, denoted by 1-1001-1 (see Eq. 11 in the manuscript)


## 1. Cα-based helix assignment rules

Helix fraction script (HelixFracDihed15.f90) is used to test different helix assignment rules

The script requires xtc library to compile:
```
ifort -o gethelixfrac HelixFracDihed15.f90 -lxdrf -L/xtc
```

#### To run the script:
**Options:**
```
-x(-xtc): xtc file
-o(-out): output file
-dihed(-dihedral): dihedral angle range
-assign: h' assignment rule
-nn(-nneighbor): number of neighboring 2*n residues
-e(-eq): equilibration step
-block: number of blocks
```
for (i,i+2) rule:
```
./gethelixfrac -x traj.xtc -o helix.dat -dihed 0.25 1.3 -assign 0110 -nn 1 -block 5 -eq 1000
```
for (i,i+4) rule:
```
./gethelixfrac -x traj.xtc -o helix.dat -dihed 0.25 1.7 -assign 0110 -nn 3 010 -block 5 -eq 1000
```
 The output is then compared with the DSSP-based helix fraction of Ala40 peptide from atomistic simulation (/2.AA-A40 directory)

## 2. Parameterization

* A40 simulations with different **eps_d** values to figure out the reference value for highest helicity to match the experimental helical propensity for Alanine. The example is given in /3.Parameterization/A40_example/

* Helix propensity code (HelixPropensity.f90) was used to compute helix propensity (w) directly from the simulation trajectory. To compile: 
```
ifort -o gethelixprop HelixPropensity.f90 -lxdrf -L/xtc
```
**Options:**
```
-x(-xtc): xtc file
-o(-out): output file
-dihed(-dihedral): dihedral angle range
-assign: h' assignment rule
-nn(-nneighbor): number of neighboring 2*n residues
-e(-eq): equilibration step
-block: number of blocks
-nmc: number of MC steps
-s(-seed): RNG seed
```
to run: 

for (i,i+2) rule:
```
./gethelixprop -x alanine.xtc -p alanine.pdb -o helicity_ala.dat -dihed 0.25 1.3 -assign 0110 -nn 1 010 -block 5 -eq 100000 -nmc 100000 -seed 1234567
```
for (i,i+4) rule:
```
./gethelixprop -x alanine.xtc -p alanine.pdb -o helicity_ala.dat -dihed 0.25 1.7 -assign 0110 -nn 3 010 -block 5 -eq 100000 -nmc 100000 -seed 1234567
```

* **eps_d** for remaining 19 residues determined using a host-guest system such as A20XA20 or A20X4A20
The example is given in parameterization/A20XA20_example/ directory

## 3. Validation (check /4.Validation directory)

* TDP-43 CTD files was added as an example for the validation step
*  Simulated helix fraction was calculated from TDP-43 CTD single chain simulation trajectory using the parameters for (i,i+4) rule
* From the available NMR data from BMRM database, a single residue-specific secondary structure propensity (SSP) score is calculated based on the deviations of NMR chemical shifts from Poulsen IDP/IUP random coil chemical shifts.
* Poulsen IDP/IUP random coil chemical shifts: https://spin.niddk.nih.gov/bax/nmrserver/Poulsen_rc_CS/
* SSP: http://pound.med.utoronto.ca/software.html

## 4. Rg calculation


## LICENSE AND DISCLAIMER

Redistribution and use of this software in source and binary forms, with or without modification, are permitted provided that this statement and the following disclaimer are retained. THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
