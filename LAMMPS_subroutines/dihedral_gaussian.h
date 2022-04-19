/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   Functional form: U(x)=-log(Ua(x)*exp(ea) + Ub(x)*exp(-ea) + Uc(x)*exp(ea-1.5) + Ud(x)*exp(-ea-1.0))
                    Ua(x)=exp(-ka*(x-fa)**2)
                    Ub(x)=exp(-kb*(x-fb)**2) + exp(-kb*(x-fb-2*pi)**2)
                    Uc(x)=exp(-kc*(x-fc)**2) + exp(-kc*(x-fc+2*pi)**2)
                    Ud(x)=exp(-kd*(x-fd)**2) + exp(-kd*(x-fd-2*pi)**2) 
                    ka = 20.0, kb = 3.0, kc = 1.0, kd = 1.0
                    fa = 0.9, fb = -1.9, fc = 0.9, fd = -3.0
------------------------------------------------------------------------- */

#ifdef DIHEDRAL_CLASS

DihedralStyle(gaussian,DihedralGaussian)

#else

#ifndef LMP_DIHEDRAL_GAUSSIAN_H
#define LMP_DIHEDRAL_GAUSSIAN_H

#include <cstdio>
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralGaussian : public Dihedral {
 public:
  DihedralGaussian(class LAMMPS *);
  ~DihedralGaussian();
  void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  double *epsdihed;

  void allocate();
};

}

#endif
#endif
