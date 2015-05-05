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

#ifndef MPS_LATTICE_H
#define MPS_LATTICE_H

#include <tensor/sparse.h>

namespace mps {

  using namespace tensor;

  /** Class representing fermionic or hard-core-bosons particles hopping in a
   * finite lattice. The lattice constructs operators representing the motion of
   * particles, their density, their interactions, assuming that there is a
   * fixed (beforehand) number of particles at all times.
   */
  class Lattice {

    const tensor::index number_of_sites;
    const int number_of_particles;
    const Indices configurations;

    static const Indices filtered_states(int sites, int number_of_particles);

  public:

    enum particle_kind_t {
      /** The lattice will contain impenetrable bosonic particles. */
      HARD_CORE_BOSONS = 0,
      /** The lattice will contain fermions (in Jordan-Wigner representation).*/
      FERMIONS = 1
    };

    /** Construct the internal representation for a lattice with N particles in
        those 'sites'*/
    Lattice(int sites, int N);

    /** Hopping operator for a particle between two sites. Returns the
        equivalent of \$a^\dagger_{to}a_{from}\$. */
    const RSparse hopping_operator(int to, int from, particle_kind_t kind = FERMIONS) const;
    /** Number operator for the given lattice site.*/
    const RSparse number_operator(int site) const;
    /** Hubbard interaction between different lattice site. It implements
        operator \$ n_{site1} n_{site2} \$.*/
    const RSparse interaction_operator(int site1, int site2) const;

    /** Full Hamiltonian containing hopping of kind and
        interactions. Matrix J(i,j) is nonzero when there is hopping between
        sites 'i' and 'j', and interactions U(i,j) among those sites
        too. Entries in these matrices can be zero. */
    const RSparse Hamiltonian(const RTensor &J, const RTensor &interactions,
			      double mu, particle_kind_t kind = FERMIONS) const;
    /** Full Hamiltonian containing hopping of kind and
        interactions. Matrix J(i,j) is nonzero when there is hopping between
        sites 'i' and 'j', and interactions U(i,j) among those sites
        too. Entries in these matrices can be zero. */
    const CSparse Hamiltonian(const CTensor &J, const CTensor &interactions,
			      double mu, particle_kind_t kind = FERMIONS) const;

    /** Number of sites in the lattice. */
    int size() const;
    /** Preconfigured number of particles. */
    int particles() const;
    /** Dimensionality of the constrained Hilbert space. */
    tensor::index dimension() const;
  };
  
}

#endif /* !MPS_LATTICE_H */