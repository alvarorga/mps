// -*- mode: c++; fill-column: 80; c-basic-offset: 2; indent-tabs-mode: nil -*-
/*
    Copyright (c) 2010 Juan Jose Garcia Ripoll

    Tensor is free software; you can redistribute it and/or modify it
    under the terms of the GNU Library General Public License as published
    by the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Library General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <mps/mps.h>
#include <mps/mps_algorithms.h>
#include <mps/mpo.h>
#include <math.h>

namespace mps {

  /*
     Computes how close is a state psi to be an eigenstate of the 
     Hamiltonian. We compute the following quantity:
         
       eig_F = <psi|H|psi>/sqrt(<psi|H^2|psi>)
     
     We will build and intermediate state: 'Hpsi' = H|psi>, which we
     simplify according to its SVD decomposition in order to reduce
     the size of the computation

     Input:
       H: Hamiltonian MPO
       psi: state whose proximity to be an eigenstate of H we want to
         measure
       ptrE: pointer to expected value <psi|H|psi>
     Simplification algorithm inputs:
       simp_tol: tolerance of the SVD singular values
       simp_Dmax: maximum bond dimension of Hpsi simplified
       simp_sweeps: number of sweeps
       simp_err: error of the simplification
     Ouput:
       eig_F: eigenstate fidelity
   */
  
  template<class MPS, class MPO>
  static double
  do_eigenstate_fidelity(const MPO &H, const MPS &psi, 
                         double &simp_err, double simp_tol,
                         tensor::index simp_sweeps, tensor::index simp_Dmax,
                         double *ptrE)
  {
    MPS Hpsi(apply(H, psi));

    // Simplify Hpsi
    int simp_sense = -1;
    bool normalize = 0;
    simp_err = simplify_obc(&Hpsi, Hpsi, &simp_sense, simp_sweeps,
                            normalize, simp_Dmax, simp_tol);

    // Compute <psi|H|psi> if it has not been initialized
    double E;
    if (ptrE)
        E = *ptrE;
    else
        E = real(scprod(psi, Hpsi));

    // Compute eigenstate fidelity
    double eig_F = tensor::abs(E)/real(sqrt(scprod(Hpsi, Hpsi)));
    return eig_F;
  }

} // namespace mps
