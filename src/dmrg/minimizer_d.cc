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

#include "minimizer.cc"

namespace mps {

  double minimize(const RMPO &H, RMPS *psi, const MinimizerOptions &opt,
                  const RMPO &constraints, double value,
                  double &eig_fidelity, double &simp_err)
  {
    Minimizer<RMPO> min(opt, H, *psi);
    min.add_constraint(constraints, value);
    return min.full_sweep(psi, eig_fidelity, simp_err);
  }

  double minimize(const RMPO &H, RMPS *psi, const MinimizerOptions &opt,
                  const RMPO &constraints, double value)
  {
    Minimizer<RMPO> min(opt, H, *psi);
    min.add_constraint(constraints, value);
    double eig_fidelity = -1.;
    double simp_err = -1.;
    return min.full_sweep(psi, eig_fidelity, simp_err);
  }

  double minimize(const RMPO &H, RMPS *psi, const MinimizerOptions &opt,
                  double &eig_fidelity, double &simp_err)
  {
    Minimizer<RMPO> min(opt, H, *psi);
    return min.full_sweep(psi, eig_fidelity, simp_err);
  }

  double minimize(const RMPO &H, RMPS *psi, const MinimizerOptions &opt)
  {
    Minimizer<RMPO> min(opt, H, *psi);
    double eig_fidelity =- 1.;
    double simp_err = -1.;
    return min.full_sweep(psi, eig_fidelity, simp_err);
  }

  double minimize(const RMPO &H, RMPS *psi,
                  double &eig_fidelity, double &simp_err)
  {
    return minimize(H, psi, MinimizerOptions(), eig_fidelity, simp_err);
  }

  double minimize(const RMPO &H, RMPS *psi)
  {
    double eig_fidelity = -1.;
    double simp_err = -1.;
    return minimize(H, psi, MinimizerOptions(), eig_fidelity, simp_err);
  }

} // namespace mps
