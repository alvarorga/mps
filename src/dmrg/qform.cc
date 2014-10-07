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

#include <algorithm>
#include <mps/qform.h>
#include <mps/mps_algorithms.h>

namespace mps {

  template<class MPO>
  QuadraticForm<MPO>::QuadraticForm(const MPO &mpo, const MPS &bra, const MPS &ket, int start) :
    matrix_(make_matrix_database(mpo)),
    pairs_(make_pairs(mpo))
  {
    if (bra[0].dimension(0) != 1 || ket[0].dimension(0) != 1) {
      std::cerr << "Due to laziness of their programmers, mps does not implement QForm for PBC";
      abort();
    }
    if (start == 0) {
      // Prepare the right matrices from site size()-1 to 0, incrementally
      // There is only one left matrix, at 0, which is empty.
      current_site_ = size() - 1;
      while (here() != 0)
	propagate_left(bra[here()], ket[here()]);
    } else {
      current_site_ = 0;
      while (here() != size()-1)
	propagate_right(bra[here()], ket[here()]);
    }
  }

  template<class MPO>
  typename QuadraticForm<MPO>::matrix_database_t
  QuadraticForm<MPO>::make_matrix_database(const MPO &mpo)
  {
    // We only support open boundary condition problems
    assert(mpo[0].dimension(0));
    index L = mpo.size();
    matrix_database_t output(L+1, matrix_array_t(1, elt_t()));
    for (index i = 0; i < L; i++) {
      output.at(i) = matrix_array_t(mpo[i].dimension(0), elt_t());
    }
    return output;
  }

  template<class MPO>
  typename QuadraticForm<MPO>::pair_tree_t
  QuadraticForm<MPO>::make_pairs(const MPO &mpo)
  {
    pair_tree_t output(mpo.size());
    for (index m = 0; m < mpo.size(); m++) {
      const elt_t &tensor = mpo[m];
      for (index i = 0; i < tensor.dimension(0); i++) {
	for (index j = 0; j < tensor.dimension(3); j++) {
	  Pair p(i, j, tensor);
	  if (!p.is_empty())
	    output.at(m).push_back(p);
	}
      }
    }
    return output;
  }

  template<class tensor>
  static void maybe_add(tensor *a, const tensor &b)
  {
    *a = (a->is_empty())? b : (*a + b);
  }

  template<class MPO>
  void QuadraticForm<MPO>::propagate_left(const elt_t &braP, const elt_t &ketP)
  {
    assert(here() >= 0);
    const matrix_array_t &mr = right_matrices(here());
    matrix_array_t &new_mr = right_matrices(here());
    std::fill(new_mr.begin(), new_mr.end(), elt_t());
    for (pair_iterator_t it = pairs_[here()].begin(), end = pairs_[here()].end();
	 it != end;
	 it++)
      {
	maybe_add<elt_t>(&new_mr.at(it->right_ndx),
			 prop_matrix(mr[it->left_ndx], -1, braP, ketP, &it->op));
      }
    --current_site_;
  }

  template<class MPO>
  void QuadraticForm<MPO>::propagate_right(const elt_t &braP, const elt_t &ketP)
  {
    assert(here()+1 < size());
    const matrix_array_t &ml = left_matrices(here());
    matrix_array_t &new_ml = left_matrices(here()+1);
    std::fill(new_ml.begin(), new_ml.end(), elt_t());
    for (pair_iterator_t it = pairs_[here()].begin(), end = pairs_[here()].end();
	 it != end;
	 it++)
      {
	maybe_add<elt_t>(&new_ml.at(it->right_ndx),
			 prop_matrix(ml[it->left_ndx], +1, braP, ketP, &it->op));
      }
    ++current_site_;
  }

  template<class elt_t>
  static elt_t compose(const elt_t &L, const elt_t &op, const elt_t &R)
  {
    index a1,a2,b1,b2,a3,b3;
    // L(a1,a2,b1,b2) op(i,j) R(a3,a1,b3,b1) -> H([a2,i,a3],[b2,j,b3])
    L.get_dimensions(&a1, &a2, &b1, &b2);
    R.get_dimensions(&a3, &a1, &b3, &b1);
    assert(a1 == 1 && b1 == 1);
    // Remember that kron(A(i,j),B(k,l)) -> C([k,i],[l,j])
    return kron(kron(reshape(R, a3,b3), op), reshape(L, a2,b2));
  }

  template<class elt_t>
  static elt_t compose(const elt_t &L, const elt_t &op1, const elt_t &op2, const elt_t &R)
  {
    // L(a1,a2,b1,b2) op(i,j) R(a3,a1,b3,b1) -> H([a2,i,a3],[b2,j,b3])
    index a1,a2,b1,b2,a3,b3;
    L.get_dimensions(&a1, &a2, &b1, &b2);
    R.get_dimensions(&a3, &a1, &b3, &b1);
    assert(a1 == 1 && b1 == 1);
    // Remember that kron(A(i,j),B(k,l)) -> C([k,i],[l,j])
    return kron(kron(kron(reshape(R, a3,b3), op2), op1), reshape(L, a2,b2));
  }

  template<class MPO>
  const typename QuadraticForm<MPO>::elt_t
  QuadraticForm<MPO>::single_site_matrix()
  {
    elt_t output;
    for (pair_iterator_t it = pairs_[here()].begin(), end = pairs_[here()].end();
	 it != end;
	 it++)
      {
	const elt_t &vl = left_matrix(here(), it->left_ndx);
	const elt_t &vr = right_matrix(here(), it->right_ndx);
	if (!vl.is_empty() && !vr.is_empty()) {
	  maybe_add<elt_t>(&output, compose(vl, it->op, vr));
	}
      }
    return output;
  }


  template<class MPO>
  const typename QuadraticForm<MPO>::elt_t
  QuadraticForm<MPO>::two_site_matrix()
  {
    assert(here() + 1 < size());
    elt_t output;
    for (pair_iterator_t it1 = pairs_[here()].begin(), end1 = pairs_[here()].end();
	 it1 != end1;
	 it1++)
      {
	for (pair_iterator_t it2 = pairs_[here()+1].begin(), end2 = pairs_[here()+1].end();
	     it2 != end2;
	     it2++)
	  if (it1->right_ndx == it2->left_ndx) {
	    const elt_t &vl = left_matrix(here(), it1->left_ndx);
	    const elt_t &vr = right_matrix(here()+1, it2->right_ndx);
	    if (!vl.is_empty() && !vr.is_empty()) {
	      maybe_add(&output, compose(vl, it1->op, it2->op, vr));
	    }
	  }
      }
    return output;
  }

}
