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

#include <mps/mpo.h>
#include "eigenstate_fidelity.cc"

namespace mps {

  double
  eigenstate_fidelity(const CMPO &H, const CMPS &psi,
                      double &simp_err, double simp_tol,
                      tensor::index simp_sweeps, tensor::index simp_Dmax,
                      double *ptrE)
  {
    return do_eigenstate_fidelity(H, psi, simp_err, simp_tol,
                                  simp_sweeps, simp_Dmax, ptrE);
  }
  
} // namespace mps
