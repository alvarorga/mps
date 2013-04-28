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

#ifndef MPDO_MPDO_H
#define MPO_MPO_H

#include <mps/mps.h>
#include <mps/rmpo.h>
#include <mps/cmpo.h>

namespace mps {

  using namespace tensor;

  void add_local_term(RMPO &mpdo, const RTensor &Hloc);

  void add_interaction(RMPO &mpdo, const RTensor &Hi, const RTensor &Hj);

  void add_local_term(CMPO &mpdo, const CTensor &Hloc);

  void add_interaction(CMPO &mpdo, const CTensor &Hi, const CTensor &Hj);

  RMPS apply(const RMPO &mpdo, const RMPS &state);

  CMPS apply(const CMPO &mpdo, const CMPS &state);

  double expected(const RMPS &bra, const RMPO &op, const RMPS &ket);

  double expected(const RMPS &bra, const RMPO &op);

  cdouble expected(const CMPS &bra, const CMPO &op, const CMPS &ket);

  cdouble expected(const CMPS &bra, const CMPO &op);

}

#endif /* !TENSOR_MPO_H */
