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
#include <mps/lform.h>
#include <mps/mps_algorithms.h>

namespace mps {

  template<class MPS>
  LinearForm<MPS>::LinearForm(const MPS &bra, const MPS &ket, int start) :
    weight_(tensor_t::ones(igen << 1)),
    bra_(std::vector<MPS>(1,bra)),
    matrix_(make_matrix_array())
  {
    initialize_matrices(start, ket);
    // for (int i = 0; i < ket.size(); i++)
    //   std::cout << "\nL[" << i << "]=" << left_matrix(0,i) << std::endl
    //             << "Q[" << i << "]=" << bra[i] << std::endl
    //             << "P[" << i << "]=" << ket[i] << std::endl
    //             << "R[" << i << "]=" << right_matrix(0,i) << std::endl;
  }

  template<class MPS>
  LinearForm<MPS>::LinearForm(const tensor_t &weight, const std::vector<MPS> &bras,
			      const MPS &ket, int start) :
    weight_(weight), bra_(bras),
    matrix_(make_matrix_array())
  {
    initialize_matrices(start, ket);
  }

  template<class MPS>
  void LinearForm<MPS>::initialize_matrices(int start, const MPS &ket)
  {
    for (current_site_ = 0; here() < start; ) {
      propagate_right(ket[here()]);
    }
    for (current_site_ = size() - 1; here() > start; ) {
      propagate_left(ket[here()]);
    }
  }

  template<class MPS>
  typename LinearForm<MPS>::matrix_database_t
  LinearForm<MPS>::make_matrix_array()
  {
    return matrix_database_t(number_of_bras(), matrix_array_t(size()+1, tensor_t()));
  }

  template<class tensor>
  static void maybe_add(tensor *a, const tensor &b)
  {
    if (a->is_empty())
      *a = b;
    else
      *a += b;
  }

  template<class MPS>
  void LinearForm<MPS>::propagate(const tensor_t &ketP, int sense)
  {
    if (sense > 0) propagate_right(ketP); else propagate_left(ketP);
  }

  template<class MPS>
  void LinearForm<MPS>::propagate_left(const tensor_t &ketP)
  {
    if (here() == 0)
      return;
    for (int i = 0; i < number_of_bras(); i++) {
      // std::cout << "PL @ " << i << std::endl
      //           << "Q=" << bra_[i][here()] << std::endl
      //           << "P=" << ketP << std::endl
      //           << "R=" << right_matrix(i,here())
      //           << std::endl;
      right_matrix(i, here()-1) =
	prop_matrix(right_matrix(i, here()), -1, bra_[i][here()], ketP);
    }
    --current_site_;
  }

  template<class MPS>
  void LinearForm<MPS>::propagate_right(const tensor_t &ketP)
  {
    if (here()+1 == size())
      return;
    for (int i = 0; i < number_of_bras(); i++) {
      // std::cout << "PR @ " << i << std::endl
      //           << "L=" << left_matrix(i,here()) << std::endl
      //           << "Q=" << bra_[i][here()] << std::endl
      //           << "P=" << ketP << std::endl;
      left_matrix(i, here()+1) =
	prop_matrix(left_matrix(i, here()), +1, bra_[i][here()], ketP);
    }
    ++current_site_;
  }

  template<class tensor_t>
  static tensor_t compose(const tensor_t &L, const tensor_t &P, const tensor_t &R)
  {
    if (L.is_empty()) {
      return compose(tensor_t::ones(1,1,1,1), P, R);
    }
    if (R.is_empty()) {
      return compose(L, P, tensor_t::ones(1,1,1,1));
    }
    // std::cout << " L=" << L << std::endl
    //           << " P=" << P << std::endl
    //           << " R=" << R << std::endl;
    index a1,a2,b1,b2,a3,b3,i;
    L.get_dimensions(&a1, &b1, &a2, &b2);
    R.get_dimensions(&a3, &b3, &a1, &b1);
    P.get_dimensions(&a2, &i, &a3);
    assert(a1 == 1 && b1 == 1);
    // Reshape L -> L(a2,b2), R -> R(a3,b3)
    // and return L(a2,b2) P(a2,i,a3) R(a3,b3)
    return fold(fold(reshape(L, a2,b2), 0, P, 0), -1,
		     reshape(R, a3, b3), 0);
  }

  template<class MPS>
  const typename LinearForm<MPS>::tensor_t
  LinearForm<MPS>::single_site_vector() const
  {
    tensor_t output;
    for (int i = 0; i < number_of_bras(); i++) {
      maybe_add(&output,
		compose<tensor_t>(left_matrix(i, here()), tensor::conj(weight_[i] * bra_[i][here()]),
                                  right_matrix(i, here())));
    }
    return output;
  }

  template<class tensor_t>
  static tensor_t compose4(const tensor_t L, const tensor_t &P1, const tensor_t &P2, const tensor_t &R)
  {
    if (L.is_empty()) {
      return compose4(tensor_t::ones(1,1,1,1), P1, P2, R);
    }
    if (R.is_empty()) {
      return compose4(L, P1, P2, tensor_t::ones(1,1,1,1));
    }
    index a1,a2,b1,b2,a3,b3,a4,b4,i,j;
    L.get_dimensions(&a1, &b1, &a2, &b2);
    R.get_dimensions(&a4, &b4, &a1, &b1);
    assert(a1 == 1 && b1 == 1);
    P1.get_dimensions(&a2, &i, &a3);
    P2.get_dimensions(&a3, &j, &a4);
    // P(a2,i,j,a4) = P1(a2,i,a3)P2(a3,j,a4)
    tensor_t P = fold(P1, -1, P2, 0);
    // Reshape L -> L(a2,b2), R -> R(a4,b4)
    // and return L(a2,b2) P(a2,i,i,a4) R(a4,b4)
    return fold(fold(reshape(L, a2,b2), 0, P, 0), -1,
		reshape(R, a4, b4), 0);
  }

  template<class MPS>
  const typename LinearForm<MPS>::tensor_t
  LinearForm<MPS>::two_site_vector(int sense) const
  {
    tensor_t output;
    index i, j;
    if (sense > 0) {
      i = here();
      j = i+1;
      assert(j < size());
    } else {
      j = here();
      assert(j > 0);
      i = j - 1;
    }
    for (int n = 0; n < number_of_bras(); n++) {
      maybe_add(&output,
		compose4<tensor_t>(left_matrix(n, i),
                                   tensor::conj(weight_[n] * bra_[n][i]),
                                   tensor::conj(bra_[n][j]),
                                   right_matrix(n, j)));
    }
    return output;
  }

  template<class MPS>
  double
  LinearForm<MPS>::norm2() const
  {
    number_t x, v = number_zero<number_t>();
    for (int i = 0; i < number_of_bras(); i++) {
      for (int j = 0; j <= i; j++) {
        x = tensor::conj(weight_[i]) * weight_[j] * scprod(bra_[i], bra_[j]);
        v += x;
        if (i != j) v += tensor::conj(x);
      }
    }
    return sqrt(tensor::abs(v));
  }

}
