/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(bond/exchange/assoc,FixBondExchangeAssoc);
// clang-format on
#else

#ifndef LMP_FIX_BOND_EXCHANGE_ASSOC_H
#define LMP_FIX_BOND_EXCHANGE_ASSOC_H

#include "fix.h"

#include <vector>

namespace LAMMPS_NS {

class FixBondExchangeAssoc : public Fix {
 public:
  FixBondExchangeAssoc(class LAMMPS *, int, char **);
  ~FixBondExchangeAssoc() override;

  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void post_integrate() override;
  int modify_param(int, char **) override;
  double compute_scalar() override;
  double compute_vector(int) override;
  double memory_usage() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  int ntry, seed;
  int tflag, topology_checked, defensive_checks;
  int pivot_type, free_type, bound_type, bond_type;
  double rmin, rmax, rminsq, rmaxsq, pkin;

  char *id_temp;
  class Compute *temperature;
  class RanMars *random;
  class NeighList *list;

  bigint ntrials, naccept, nreject, n_skip_pkin;
  bigint n_no_pivot, n_no_leaving, n_no_attacker;

  int *pivot_counts;
  int ntotal_pivots;

  std::vector<int> pivots;
  std::vector<int> attackers;
  std::vector<int> active_ranks;
  std::vector<tagint> topology_send;
  std::vector<tagint> topology_recv;

  void build_pivots();
  int choose_active_rank();
  int pick_random(const std::vector<int> &);

  bool is_pivot(int) const;
  bool is_free_A(int) const;
  bool is_bound_A(int) const;
  bool in_group(int) const;
  bool in_window(int, int) const;
  double dist_rsq(int, int) const;
  bool stores_bond(int, tagint, int) const;
  bool has_bond(int, int, int) const;

  void check_pair_symmetry();
  void check_lj_cut_pair_symmetry();
  void check_sampled_pair_symmetry();
  void check_initial_topology();

  int find_bound_partner(int);
  void build_attackers(int, std::vector<int> &);

  double bond_eng(int, int, int) const;
  double pair_eng(int, int, int, int, int) const;
  double compute_bond_delta_energy(int, int, int) const;
  double compute_pair_delta_energy(int, int, int) const;
  bool accept_move(double, double);

  void commit_move_by_tags(tagint, tagint, tagint);
  void sync_topology_all();
  void sync_topology_for_move(tagint, tagint, tagint);
  int topology_record_size() const;
  void pack_owned_topology_record(tagint, tagint *) const;
  void unpack_topology_record(tagint, const tagint *);
};

}    // namespace LAMMPS_NS

#endif
#endif
