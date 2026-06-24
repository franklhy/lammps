// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   fix bond/exchange/assoc: rank-local associative dynamic bond exchange

   This first implementation performs one globally serialized rank-local MC
   stream.  A trial is attempted on one MPI rank at a time.  The pivot B is
   owned by the active rank; the leaving and attacking A atoms may be owned or
   ghost atoms in that rank's local view:

     Ai-B + Aj -> Ai + B-Aj

   The old and new B-A distances must both be inside the same reactive
   window, which makes the reverse proposal immediately accessible.  The
   simple symmetric version assumes both A and B are monovalent with respect
   to the exchanged dynamic bond.  This makes the forward and reverse proposal
   counts equal, so the proposal ratio is unity.  The energy model is
   intentionally restricted to short-range two-body pair styles with
   Pair::single(), bond styles with Bond::single(), no kspace, and no
   angles/dihedrals/impropers.  1-3 and 1-4 special-bond weights must be unity
   so the only topology-dependent pair correction is the exchanged 1-2
   relation.
------------------------------------------------------------------------- */

#include "fix_bond_exchange_assoc.h"

#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "random_mars.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Parse the fix command, validate user arguments, initialize output flags,
   and create a rank-local random stream.  Physics restrictions that require
   fully initialized force styles are checked later in init().
------------------------------------------------------------------------- */

FixBondExchangeAssoc::FixBondExchangeAssoc(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_temp(nullptr), temperature(nullptr), random(nullptr),
  list(nullptr), pivot_counts(nullptr), ntotal_pivots(0)
{
  std::string fixname = fmt::format("fix {}", style);
  if (narg < 18) utils::missing_cmd_args(FLERR, fixname, error);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  ntry = utils::inumeric(FLERR, arg[4], false, lmp);
  seed = utils::inumeric(FLERR, arg[5], false, lmp);
  rmin = utils::numeric(FLERR, arg[6], false, lmp);
  rmax = utils::numeric(FLERR, arg[7], false, lmp);

  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix {} Nevery value {}", style, nevery);
  if (ntry <= 0) error->all(FLERR, 4, "Illegal fix {} Ntry value {}", style, ntry);
  if (seed <= 0) error->all(FLERR, 5, "Illegal fix {} seed value {}", style, seed);
  if (rmin < 0.0 || rmax < rmin)
    error->all(FLERR, "Illegal fix {} reactive window [{},{}]", style, rmin, rmax);

  rminsq = rmin * rmin;
  rmaxsq = rmax * rmax;

  pivot_type = free_type = bound_type = bond_type = 0;
  pkin = -1.0;
  tflag = 0;
  topology_checked = 0;
  defensive_checks = 0;

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "pivot_type") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, fixname + " pivot_type", error);
      pivot_type = utils::expand_type_int(FLERR, arg[iarg + 1], Atom::ATOM, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "free_type") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, fixname + " free_type", error);
      free_type = utils::expand_type_int(FLERR, arg[iarg + 1], Atom::ATOM, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bound_type") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, fixname + " bound_type", error);
      bound_type = utils::expand_type_int(FLERR, arg[iarg + 1], Atom::ATOM, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bond_type") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, fixname + " bond_type", error);
      bond_type = utils::expand_type_int(FLERR, arg[iarg + 1], Atom::BOND, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pkin") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, fixname + " pkin", error);
      pkin = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "defensive_checks") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, fixname + " defensive_checks", error);
      defensive_checks = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else error->all(FLERR, iarg, "Unknown fix {} keyword {}", style, arg[iarg]);
  }

  if (pivot_type < 1 || pivot_type > atom->ntypes)
    error->all(FLERR, "Invalid pivot_type in fix {} command", style);
  if (free_type < 1 || free_type > atom->ntypes)
    error->all(FLERR, "Invalid free_type in fix {} command", style);
  if (bound_type < 1 || bound_type > atom->ntypes)
    error->all(FLERR, "Invalid bound_type in fix {} command", style);
  if (pivot_type == free_type || pivot_type == bound_type || free_type == bound_type)
    error->all(FLERR, "Fix {} requires distinct pivot_type, free_type, and bound_type labels",
               style);
  if (bond_type < 1 || bond_type > atom->nbondtypes)
    error->all(FLERR, "Invalid bond_type in fix {} command", style);
  if (pkin <= 0.0 || pkin > 1.0)
    error->all(FLERR, "Fix {} requires 0 < pkin <= 1", style);

  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR, Error::NOLASTLINE, "Cannot use fix {} with non-molecular systems", style);

  dynamic_group_allow = 1;
  force_reneighbor = 1;
  next_reneighbor = -1;

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 8;
  global_freq = 1;
  extscalar = 0;
  extvector = 0;

  random = new RanMars(lmp, seed + comm->me);

  // Create the default temperature compute used for the Metropolis test.
  // Users can replace it with "fix_modify ID temp compute-ID".

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp", id_temp));
  tflag = 1;

  ntrials = naccept = nreject = n_skip_pkin = 0;
  n_no_pivot = n_no_leaving = n_no_attacker = 0;

  memory->create(pivot_counts, comm->nprocs, "bond/exchange/assoc:pivot_counts");
}

/* ----------------------------------------------------------------------
   Free fix-owned resources.
------------------------------------------------------------------------- */

FixBondExchangeAssoc::~FixBondExchangeAssoc()
{
  delete random;
  if (tflag) modify->delete_compute(id_temp);
  delete[] id_temp;
  memory->destroy(pivot_counts);
}

/* ----------------------------------------------------------------------
   Request invocation after MD integration on timesteps divisible by Nevery.
   Topology changes are made before the next neighbor rebuild.
------------------------------------------------------------------------- */

int FixBondExchangeAssoc::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ----------------------------------------------------------------------
   Validate the temperature compute and force-field support.  This
   implementation assumes short-range two-body pair styles with Pair::single(),
   no kspace, no many-body terms, and no angle/dihedral/improper topology
   changes.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::init()
{
  temperature = modify->get_compute_by_id(id_temp);
  if (!temperature) {
    error->all(FLERR, "Temperature compute ID {} for fix {} does not exist", id_temp, style);
  } else {
    if (temperature->tempflag == 0)
      error->all(FLERR, "Compute ID {} for fix {} does not compute a temperature", id_temp,
                 style);
  }

  if (atom->tag_enable == 0)
    error->all(FLERR, "Fix {} requires atoms to have IDs", style);
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix {} requires an atom map", style);
  if (atom->num_bond == nullptr || atom->bond_atom == nullptr || atom->bond_type == nullptr)
    error->all(FLERR, "Fix {} requires an atom style with bonds", style);
  if (force->pair == nullptr || force->bond == nullptr)
    error->all(FLERR, "Fix {} requires pair and bond styles", style);
  if (force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support fix {}", style);
  if (force->pair->manybody_flag)
    error->all(FLERR, "Fix {} does not support many-body pair styles", style);
  if (force->kspace)
    error->all(FLERR, "Fix {} does not support kspace/long-range electrostatics", style);
  if (atom->nangles || atom->ndihedrals || atom->nimpropers)
    error->all(FLERR, "Fix {} does not yet support systems with angles, dihedrals, or impropers",
               style);
  if (force->special_lj[2] != 1.0 || force->special_coul[2] != 1.0 ||
      force->special_lj[3] != 1.0 || force->special_coul[3] != 1.0)
    error->all(FLERR, "Fix {} requires special_bonds 1-3 and 1-4 weights to be 1.0", style);
  if (atom->rmass == nullptr && atom->mass &&
      atom->mass[free_type] != atom->mass[bound_type])
    error->all(FLERR, "Fix {} requires equal masses for free_type and bound_type", style);

  comm_forward = topology_record_size();

  // A full occasional neighbor list with cutoff Rmax lets attacker search visit
  // only atoms near the chosen pivot instead of scanning all local+ghost atoms.
  auto *req = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
  req->set_cutoff(rmax);

  check_pair_symmetry();
}

/* ----------------------------------------------------------------------
   Receive the occasional full neighbor list used for attacking-A searches.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   Main MC driver for the associative exchange move:

     Ai(bound)-B + Aj(free) -> Ai(free) + B-Aj(bound)

   Algorithm, for every Nevery timestep:

   1. Refresh ghost coordinates, then use fix-specific communication to copy
      owned atom type/bond/special topology onto ghost copies.  Normal
      forward_comm() updates positions; it does not refresh molecular
      topology for ghosts.
   2. Build the owned pivot list and gather pivot counts once for this
      invocation, then build an occasional full neighbor list for local
      attacking-A searches.  Pivot labels and coordinates do not change during
      the MC sequence.
   3. Rank 0 builds and broadcasts the Ntry-trial schedule: for each trial, the
      entry is either a pkin-skip marker or the selected active rank.  If the
      trial fails this random gate, no proposal is generated.  This is
      equivalent to multiplying the Metropolis probability by pkin.
   4. If the pkin gate passes, the active MPI rank was chosen by drawing
      uniformly from the global pivot list.  Only that rank generates a
      proposal, which avoids conflicting topology edits in this rank-local
      implementation.
   5. The active rank chooses:
      a. an owned pivot B of pivot_type in the fix group,
      b. the leaving A atom i of bound_type bonded to B by bond_type,
      c. the attacking A atom j of free_type that is not directly bonded to a
         pivot.
      The leaving and attacking atoms may be owned or ghost atoms in the
      active rank's local view.  Both Ai-B and B-Aj distances must lie in the
      same [Rmin,Rmax] reactive window so the reverse move is accessible.
   6. The active rank computes dU locally.  A startup check requires free_type
      and bound_type to have identical pair interactions with every atom type,
      so pair dU only contains the exchanged Ai-B and B-Aj 1-2 special-bonds
      corrections.
   7. The active rank applies the Metropolis criterion using the current
      temperature from the fix temperature compute.  The simple symmetric
      monovalent proposal has qrev/qfwd = 1.
   8. The active rank broadcasts whether the trial was accepted.  If accepted,
      it also broadcasts the accepted atom tags, and all ranks call the
      tag-based commit routine.  Each rank updates only atoms it owns, then
      the owning ranks explicitly forward the changed type/bond/special
      topology to ghost copies before the next MC trial.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::post_integrate()
{
  if (update->ntimestep % nevery) return;

  // Ensure the active rank can use current ghost positions.  The following
  // fix-specific sync refreshes ghost type/topology, which ordinary atom
  // forward communication does not provide for molecular atom styles.

  comm->forward_comm();
  sync_topology_all();

  if (!topology_checked) {
    check_initial_topology();
    topology_checked = 1;
  }
  build_pivots();
  neighbor->build_one(list);

  double t_current = temperature->compute_scalar();
  if (t_current <= 0.0)
    error->all(FLERR, "Current temperature for fix {} is not positive", style);

  // The schedule can be precomputed because this fix changes A state labels and
  // B-A bonds, but it does not change which atoms are pivot_type or which ranks
  // own them during this invocation.  Thus the pivot counts gathered by
  // build_pivots() remain valid for all trials, and rank-weighted scheduling is
  // equivalent to drawing a globally uniform pivot one trial at a time.
  //
  // Each entry is the MPI rank that proposes that trial.  Negative entries are
  // skipped trials: -1 means no pivot exists, -2 means the pkin gate skipped it.
  active_ranks.resize(ntry);
  if (comm->me == 0) {
    for (int itrial = 0; itrial < ntry; itrial++) {
      ntrials++;
      if (random->uniform() > pkin) {
        n_skip_pkin++;
        active_ranks[itrial] = -2;
      } else {
        active_ranks[itrial] = choose_active_rank();
        if (active_ranks[itrial] < 0) n_no_pivot++;
      }
    }
  }
  MPI_Bcast(active_ranks.data(), ntry, MPI_INT, 0, world);

  int accepted_global = 0;

  for (int itrial = 0; itrial < ntry; itrial++) {
    int active_rank = active_ranks[itrial];
    if (active_rank < 0) continue;

    int accepted = 0;
    int proposed = 0;
    tagint proposal[3] = {0, 0, 0};
    double energy_delta = 0.0;

    // The active rank performs the full trial locally.  Other ranks only need
    // to know whether an accepted topology change must be committed.

    if (comm->me == active_rank) {
      int pivot = pick_random(pivots);
      int leave = find_bound_partner(pivot);

      if (leave < 0) {
        n_no_leaving++;
      } else {
        build_attackers(pivot, attackers);

        if (attackers.empty()) {
          n_no_attacker++;
        } else {
          int attack = pick_random(attackers);
          double bond_delta = compute_bond_delta_energy(pivot, leave, attack);
          double pair_delta = compute_pair_delta_energy(pivot, leave, attack);
          energy_delta = bond_delta + pair_delta;
          proposal[0] = atom->tag[pivot];
          proposal[1] = atom->tag[leave];
          proposal[2] = atom->tag[attack];
          proposed = 1;
        }
      }

      if (proposed) {
        // only the proposing rank draws the accept/reject random number

        if (accept_move(energy_delta, t_current)) {
          naccept++;
          accepted = 1;
        } else {
          nreject++;
        }
      }
    }

    MPI_Bcast(&accepted, 1, MPI_INT, active_rank, world);
    if (accepted) {
      accepted_global = 1;

      // accepted move tags are the only proposal data needed by non-active ranks

      MPI_Bcast(proposal, 3, MPI_LMP_TAGINT, active_rank, world);

      // tag-based commit lets the owners of B, i, and j update their topology

      commit_move_by_tags(proposal[0], proposal[1], proposal[2]);
      sync_topology_for_move(proposal[0], proposal[1], proposal[2]);
    }
  }

  if (accepted_global) next_reneighbor = update->ntimestep;
}

/* ----------------------------------------------------------------------
   Build the list of eligible pivot atoms owned by this rank.  The pivot is
   the only atom in a rank-local proposal that must be owned by the active
   rank; leaving and attacking atoms can be ghosts in that rank's view.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::build_pivots()
{
  pivots.clear();

  for (int i = 0; i < atom->nlocal; i++)
    if (is_pivot(i)) pivots.push_back(i);

  int nlocal_pivots = static_cast<int>(pivots.size());
  MPI_Allgather(&nlocal_pivots, 1, MPI_INT, pivot_counts, 1, MPI_INT, world);

  ntotal_pivots = 0;
  for (int i = 0; i < comm->nprocs; i++)
    ntotal_pivots += pivot_counts[i];
}

/* ----------------------------------------------------------------------
   Select the rank that owns the next pivot by drawing uniformly from the global
   list of owned pivots.  This is called on rank 0 after build_pivots() has
   gathered pivot counts.  It avoids biasing pivots on sparsely populated ranks.
   Returns -1 if no eligible pivot exists anywhere.
------------------------------------------------------------------------- */

int FixBondExchangeAssoc::choose_active_rank()
{
  if (comm->me != 0)
    error->one(FLERR, "Internal error: fix {} choose_active_rank called off rank 0", style);

  int active = -1;
  if (ntotal_pivots > 0) {
    int target = static_cast<int>(random->uniform() * ntotal_pivots);
    int offset = 0;
    for (int i = 0; i < comm->nprocs; i++) {
      if (target < offset + pivot_counts[i]) {
        active = i;
        break;
      }
      offset += pivot_counts[i];
    }
  }
  return active;
}

/* ----------------------------------------------------------------------
   Return a uniformly selected entry from a non-empty candidate vector.
------------------------------------------------------------------------- */

int FixBondExchangeAssoc::pick_random(const std::vector<int> &list)
{
  int n = static_cast<int>(list.size());
  return list[static_cast<int>(random->uniform() * n)];
}

/* ----------------------------------------------------------------------
   True if atom i is an eligible owned/ghost pivot candidate by type and group.
   build_pivots() applies the additional ownership restriction.
------------------------------------------------------------------------- */

bool FixBondExchangeAssoc::is_pivot(int i) const
{
  return in_group(i) && atom->type[i] == pivot_type;
}

/* ----------------------------------------------------------------------
   True if atom i has the free A state label and belongs to the fix group.
------------------------------------------------------------------------- */

bool FixBondExchangeAssoc::is_free_A(int i) const
{
  return in_group(i) && atom->type[i] == free_type;
}

/* ----------------------------------------------------------------------
   True if atom i has the bound A state label.  The helper is type-only;
   find_bound_partner() applies the group check required for the symmetric
   monovalent proposal.
------------------------------------------------------------------------- */

bool FixBondExchangeAssoc::is_bound_A(int i) const
{
  return atom->type[i] == bound_type;
}

/* ----------------------------------------------------------------------
   True if atom i belongs to the group the fix is applied to.
------------------------------------------------------------------------- */

bool FixBondExchangeAssoc::in_group(int i) const
{
  return atom->mask[i] & groupbit;
}

/* ----------------------------------------------------------------------
   Test the symmetric reactive distance window.  Both the old Ai-B and new B-Aj
   distances must pass this test to preserve immediate reverse accessibility.
------------------------------------------------------------------------- */

bool FixBondExchangeAssoc::in_window(int i, int j) const
{
  double rsq = dist_rsq(i, j);
  return rsq >= rminsq && rsq <= rmaxsq;
}

/* ----------------------------------------------------------------------
   Minimum-image squared distance between two local/ghost atoms.
------------------------------------------------------------------------- */

double FixBondExchangeAssoc::dist_rsq(int i, int j) const
{
  double delx = atom->x[i][0] - atom->x[j][0];
  double dely = atom->x[i][1] - atom->x[j][1];
  double delz = atom->x[i][2] - atom->x[j][2];
  domain->minimum_image(FLERR, delx, dely, delz);
  return delx * delx + dely * dely + delz * delz;
}

/* ----------------------------------------------------------------------
   Check whether one local/ghost atom actually stores a bond entry to partner.
   This tests the atom's own bond list only; it makes no assumption about
   newton_bond storage order.
------------------------------------------------------------------------- */

bool FixBondExchangeAssoc::stores_bond(int i, tagint partner, int btype) const
{
  for (int m = 0; m < atom->num_bond[i]; m++)
    if (atom->bond_atom[i][m] == partner && atom->bond_type[i][m] == btype) return true;

  return false;
}

/* ----------------------------------------------------------------------
   Check whether atoms i and j have a stored bond of type btype on either
   atom's local bond list.  This supports data-file bonds stored on either
   endpoint with newton_bond on, and the two-entry storage used with
   newton_bond off, as long as ghost topology has been synchronized.
------------------------------------------------------------------------- */

bool FixBondExchangeAssoc::has_bond(int i, int j, int btype) const
{
  return stores_bond(i, atom->tag[j], btype) || stores_bond(j, atom->tag[i], btype);
}

/* ----------------------------------------------------------------------
   Enforce the simple symmetric state-label model.  The only allowed pair
   energy change is the B-A 1-2 special-bonds correction caused by moving the
   dynamic bond.  Therefore free_type and bound_type must have identical pair
   interactions with every atom type.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::check_pair_symmetry()
{
  std::string pair_style = utils::strip_style_suffix(force->pair_style, lmp);
  if (pair_style == "lj/cut") check_lj_cut_pair_symmetry();
  else check_sampled_pair_symmetry();
}

/* ----------------------------------------------------------------------
   Exact coefficient check for pair_style lj/cut.  PairLJCut exposes epsilon
   and sigma through Pair::extract(); cutoffs are available through Pair::cutsq.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::check_lj_cut_pair_symmetry()
{
  int dim;
  auto **epsilon = (double **) force->pair->extract("epsilon", dim);
  if (epsilon == nullptr || dim != 2)
    error->all(FLERR, "Fix {} could not extract lj/cut epsilon coefficients", style);

  auto **sigma = (double **) force->pair->extract("sigma", dim);
  if (sigma == nullptr || dim != 2)
    error->all(FLERR, "Fix {} could not extract lj/cut sigma coefficients", style);

  auto same = [](double a, double b) {
    const double scale = MAX(1.0, MAX(fabs(a), fabs(b)));
    return fabs(a - b) <= 1.0e-12 * scale;
  };

  auto eps = [&](int i, int j) {
    if (i > j) std::swap(i, j);
    return epsilon[i][j];
  };

  auto sig = [&](int i, int j) {
    if (i > j) std::swap(i, j);
    return sigma[i][j];
  };

  auto cut2 = [&](int i, int j) {
    if (i > j) std::swap(i, j);
    return force->pair->cutsq[i][j];
  };

  for (int ctype = 1; ctype <= atom->ntypes; ctype++) {
    if (!same(eps(free_type, ctype), eps(bound_type, ctype)) ||
        !same(sig(free_type, ctype), sig(bound_type, ctype)) ||
        !same(cut2(free_type, ctype), cut2(bound_type, ctype)))
      error->all(FLERR,
                 "Fix {} requires identical lj/cut coefficients for free_type and bound_type "
                 "against atom type {}",
                 style, ctype);
  }
}

/* ----------------------------------------------------------------------
   Generic numerical fallback for pair styles without a coefficient API.  It
   compares Pair::single() for (free_type,C) and (bound_type,C), and for the
   reversed order, at Rmin, midpoint, and Rmax.  This is a practical check for
   simple two-body styles, not a symbolic proof of coefficient equality.  Only
   one rank performs the sampling, but it must be a rank with a valid atom index:
   some Pair::single() implementations ignore i/j, while others dereference
   atom arrays such as charge or radius.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::check_sampled_pair_symmetry()
{
  if (rmin <= 0.0)
    error->all(FLERR,
               "Fix {} requires Rmin > 0 for numerical pair-symmetry checks with pair style {}",
               style, force->pair_style);

  int nall = atom->nlocal + atom->nghost;
  int sample_rank_candidate = nall > 0 ? comm->me : comm->nprocs;
  int sample_rank = comm->nprocs;
  MPI_Allreduce(&sample_rank_candidate, &sample_rank, 1, MPI_INT, MPI_MIN, world);
  if (sample_rank == comm->nprocs)
    error->all(FLERR, "Fix {} could not sample pair interactions because no atoms exist", style);

  double rtest[3] = {rmin, 0.5 * (rmin + rmax), rmax};
  int mismatch = 0;
  const double tol = 1.0e-10;

  if (comm->me == sample_rank) {
    int isample = 0;
    for (int ctype = 1; ctype <= atom->ntypes; ctype++) {
      for (int m = 0; m < 3; m++) {
        double rsq = rtest[m] * rtest[m];
        double fpair;

        double efc = 0.0;
        if (rsq < force->pair->cutsq[free_type][ctype])
          efc = force->pair->single(isample, isample, free_type, ctype, rsq, 1.0, 1.0, fpair);
        double ebc = 0.0;
        if (rsq < force->pair->cutsq[bound_type][ctype])
          ebc = force->pair->single(isample, isample, bound_type, ctype, rsq, 1.0, 1.0, fpair);

        double scale = 1.0;
        if (fabs(efc) > scale) scale = fabs(efc);
        if (fabs(ebc) > scale) scale = fabs(ebc);
        if (fabs(efc - ebc) > tol * scale) mismatch = ctype * 10 + m + 1;

        double ecf = 0.0;
        if (rsq < force->pair->cutsq[ctype][free_type])
          ecf = force->pair->single(isample, isample, ctype, free_type, rsq, 1.0, 1.0, fpair);
        double ecb = 0.0;
        if (rsq < force->pair->cutsq[ctype][bound_type])
          ecb = force->pair->single(isample, isample, ctype, bound_type, rsq, 1.0, 1.0, fpair);

        scale = 1.0;
        if (fabs(ecf) > scale) scale = fabs(ecf);
        if (fabs(ecb) > scale) scale = fabs(ecb);
        if (fabs(ecf - ecb) > tol * scale) mismatch = ctype * 10 + m + 1;
      }
    }
  }

  int mismatch_any = 0;
  MPI_Allreduce(&mismatch, &mismatch_any, 1, MPI_INT, MPI_MAX, world);
  if (mismatch_any)
    error->all(FLERR,
               "Fix {} requires free_type and bound_type to have identical sampled pair "
               "interactions with every atom type",
               style);
}

/* ----------------------------------------------------------------------
   One-time sanity check for the monovalent state-label model.  It verifies
   the invariant used by the proposal filters: free_type atoms have no direct
   A-B bond to a pivot, bound_type atoms have exactly one dynamic bond to a
   pivot, and pivot atoms have at most one dynamic bond to a bound A.  The check
   walks 1-2 special neighbors so it does not depend on which atom stores the
   bond when newton_bond is on.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::check_initial_topology()
{
  int nall = atom->nlocal + atom->nghost;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->type[i] != free_type && atom->type[i] != bound_type &&
        atom->type[i] != pivot_type)
      continue;

    int ndynamic = 0;
    for (int m = 0; m < atom->nspecial[i][0]; m++) {
      int j = atom->map(atom->special[i][m]);
      if (j < 0 || j >= nall)
        error->one(FLERR, "Fix {} could not map a 1-2 special neighbor during topology check",
                   style);

      if ((atom->type[i] == free_type && atom->type[j] == pivot_type) ||
          (atom->type[i] == pivot_type && atom->type[j] == free_type))
        error->one(FLERR, "Fix {} found a free_type atom directly bonded to a pivot_type atom",
                   style);

      if (!has_bond(i, j, bond_type)) continue;

      if (atom->type[i] == free_type)
        error->one(FLERR, "Fix {} found a free_type atom with a dynamic B-A bond", style);

      if (atom->type[i] == bound_type) {
        if (atom->type[j] != pivot_type)
          error->one(FLERR, "Fix {} found bound_type bonded by bond_type to non-pivot atom",
                     style);
        ndynamic++;
      } else if (atom->type[i] == pivot_type) {
        if (atom->type[j] != bound_type)
          error->one(FLERR, "Fix {} found pivot_type bonded by bond_type to non-bound atom",
                     style);
        ndynamic++;
      }
    }

    if (atom->type[i] == bound_type && ndynamic != 1)
      error->one(FLERR, "Fix {} requires each bound_type atom to have one dynamic B-A bond",
                 style);
    if (atom->type[i] == pivot_type && ndynamic > 1)
      error->one(FLERR, "Fix {} requires each pivot_type atom to have at most one dynamic B-A bond",
                 style);
  }
}

/* ----------------------------------------------------------------------
   Find the leaving A partner from the pivot's 1-2 special-neighbor list.
   The special list contains bonded neighbors independent of which atom stores
   the bond when newton_bond is on.  We then verify bond_type with has_bond()
   because the special list stores connectivity but not bond types.
------------------------------------------------------------------------- */

int FixBondExchangeAssoc::find_bound_partner(int pivot)
{
  for (int m = 0; m < atom->nspecial[pivot][0]; m++) {
    int i = atom->map(atom->special[pivot][m]);
    if (i < 0 || i >= atom->nlocal + atom->nghost) continue;
    if (!is_bound_A(i)) continue;
    if (!has_bond(pivot, i, bond_type)) continue;
    if (!in_group(i)) continue;
    if (!in_window(pivot, i)) continue;
    return i;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   Build all attacking A candidates for one monovalent pivot from the occasional
   full neighbor list.  Candidates can be owned or ghost atoms, must be
   free_type and in the fix group, and must lie inside the reactive window.  The
   one-time topology check guarantees that a free_type atom has no direct A-B
   bond to a pivot.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::build_attackers(int pivot, std::vector<int> &candidates)
{
  candidates.clear();

  int *jlist = list->firstneigh[pivot];
  int jnum = list->numneigh[pivot];
  for (int jj = 0; jj < jnum; jj++) {
    int i = jlist[jj] & NEIGHMASK;
    if (!is_free_A(i)) continue;
    if (!in_window(pivot, i)) continue;
    candidates.push_back(i);
  }
}

/* ----------------------------------------------------------------------
   Bond::single() wrapper for the dynamic B-A bond energy at the current
   minimum-image distance.
------------------------------------------------------------------------- */

double FixBondExchangeAssoc::bond_eng(int btype, int i, int j) const
{
  double fbond;
  return force->bond->single(btype, dist_rsq(i, j), i, j, fbond);
}

/* ----------------------------------------------------------------------
   Pair::single() wrapper for one local/ghost pair with a prescribed special
   level: 0 = normal pair, 1 = bonded 1-2 pair.  The fix rejects non-unity 1-3
   and 1-4 weights, so no other levels are needed for the MC dU.
------------------------------------------------------------------------- */

double FixBondExchangeAssoc::pair_eng(int i, int j, int itype, int jtype,
                                      int special_level) const
{
  double rsq = dist_rsq(i, j);
  if (rsq >= force->pair->cutsq[itype][jtype]) return 0.0;

  double fpair;
  return force->pair->single(i, j, itype, jtype, rsq, force->special_coul[special_level],
                             force->special_lj[special_level], fpair);
}

/* ----------------------------------------------------------------------
   Bond part of dU for transferring the dynamic B-A bond from leave to attack.
   The old B-leave bond is removed and the new B-attack bond is created with
   the same bond type.
------------------------------------------------------------------------- */

double FixBondExchangeAssoc::compute_bond_delta_energy(int pivot, int leave, int attack) const
{
  return bond_eng(bond_type, pivot, attack) - bond_eng(bond_type, pivot, leave);
}

/* ----------------------------------------------------------------------
   Pair part of dU for the simple symmetric state-label model.  Since init()
   verifies that free_type and bound_type have identical pair interactions with
   every atom type, all type-change pair terms cancel except the changed 1-2
   special-bonds status of B-leave and B-attack.
------------------------------------------------------------------------- */

double FixBondExchangeAssoc::compute_pair_delta_energy(int pivot, int leave, int attack) const
{
  double before = 0.0;
  before += pair_eng(pivot, leave, pivot_type, bound_type, 1);
  before += pair_eng(pivot, attack, pivot_type, free_type, 0);

  double after = 0.0;
  after += pair_eng(pivot, leave, pivot_type, free_type, 0);
  after += pair_eng(pivot, attack, pivot_type, bound_type, 1);

  return after - before;
}

/* ----------------------------------------------------------------------
   Apply the Metropolis criterion with the current temperature.  The kinetic
   prefactor pkin has already been applied as a pre-proposal random gate.
   The simple symmetric monovalent proposal has qrev/qfwd = 1.
------------------------------------------------------------------------- */

bool FixBondExchangeAssoc::accept_move(double dU, double t_current)
{
  double boltz = 1.0;
  if (dU > 0.0) boltz *= exp(-dU / (force->boltz * t_current));
  double prob = MIN(1.0, boltz);

  return random->uniform() < prob;
}

/* ----------------------------------------------------------------------
   Commit an accepted exchange by global atom tags.  Each rank modifies only
   its owned copies of B, leaving A, and attacking A.  Ghost copies are then
   refreshed explicitly by sync_topology_for_move() before another MC trial is
   attempted.  The old bond is deleted from whichever endpoint actually stores
   it, so data-file bonds may have listed either atom first.  Newly created
   bonds follow fix bond/create's convention: with newton_bond on the lower-tag
   atom stores the bond, and with newton_bond off both atoms store it.  Each
   owned endpoint also updates its 1-2 special-neighbor entry here, because the
   special list is the topology data used by pair special_bonds scaling.  If an
   owned atom stores both old and new bonds, replacing the partner in-place
   follows fix bond/swap and preserves bond-history slot bookkeeping.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::commit_move_by_tags(tagint tag_pivot, tagint tag_leave,
                                               tagint tag_attack)
{
  int pivot = atom->map(tag_pivot);
  int leave = atom->map(tag_leave);
  int attack = atom->map(tag_attack);

  // Detect where the old bond is actually stored on owned atoms.  Initial
  // data-file bonds with newton_bond on are stored on the first atom listed in
  // the Bonds section, so the old entry may be on pivot, leave, or both.  New
  // bonds created by this fix use the deterministic lower-tag convention used
  // by fix bond/create; with newton_bond off, both endpoints store the bond.

  int newton_bond = force->newton_bond;
  int old_on_pivot =
      pivot >= 0 && pivot < atom->nlocal && stores_bond(pivot, tag_leave, bond_type);
  int old_on_leave =
      leave >= 0 && leave < atom->nlocal && stores_bond(leave, tag_pivot, bond_type);
  int new_on_pivot = !newton_bond || tag_pivot < tag_attack;
  int new_on_attack = !newton_bond || tag_attack < tag_pivot;

  // Optional defensive consistency check.  A valid LAMMPS topology should have
  // one stored old-bond entry with newton_bond on and two with newton_bond off.
  // When enabled, this catches inconsistent input or bugs before the fix
  // deletes/creates bond entries.

  if (defensive_checks) {
    int old_entries_local = old_on_pivot + old_on_leave;
    int old_entries = 0;
    MPI_Allreduce(&old_entries_local, &old_entries, 1, MPI_INT, MPI_SUM, world);

    int expected_old_entries = newton_bond ? 1 : 2;
    if (old_entries != expected_old_entries)
      error->all(FLERR,
                 "Fix {} expected {} stored old dynamic bond entries but found {} for atoms {} and {}",
                 style, expected_old_entries, old_entries, tag_pivot, tag_leave);
  }

  // Check capacity only for atoms that gain a stored bond slot.  If the pivot
  // stores both old and new bonds, the old partner is replaced in-place and no
  // extra slot is needed.

  if (pivot >= 0 && pivot < atom->nlocal && !old_on_pivot && new_on_pivot &&
      atom->num_bond[pivot] >= atom->bond_per_atom)
    error->one(FLERR, "New bond from fix {} exceeded bonds per atom limit", style);
  if (attack >= 0 && attack < atom->nlocal && new_on_attack &&
      atom->num_bond[attack] >= atom->bond_per_atom)
    error->one(FLERR, "New bond from fix {} exceeded bonds per atom limit", style);

  auto histories = modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  // Delete by compacting the bond list and shifting matching bond-history slots.

  auto delete_one = [&](int i, tagint partner) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] != partner || atom->bond_type[i][m] != bond_type) continue;
      // Delete the bond between i and partner by shifting later bond slots left.
      for (int n = m; n < atom->num_bond[i] - 1; n++) {
        atom->bond_atom[i][n] = atom->bond_atom[i][n + 1];
        atom->bond_type[i][n] = atom->bond_type[i][n + 1];
        if (n_histories > 0)
          for (auto &hist : histories)
            dynamic_cast<FixBondHistory *>(hist)->shift_history(i, n, n + 1);
      }
      if (n_histories > 0)
        for (auto &hist : histories)
          dynamic_cast<FixBondHistory *>(hist)->delete_history(i, atom->num_bond[i] - 1);
      atom->num_bond[i]--;
      return true;
    }
    return false;
  };

  auto replace_one = [&](int i, tagint old_partner, tagint new_partner) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] != old_partner || atom->bond_type[i][m] != bond_type) continue;
      if (n_histories > 0)
        for (auto &hist : histories)
          dynamic_cast<FixBondHistory *>(hist)->delete_history(i, m);
      atom->bond_atom[i][m] = new_partner;
      atom->bond_type[i][m] = bond_type;
      return true;
    }
    return false;
  };

  auto add_one = [&](int i, tagint partner) {
    int m = atom->num_bond[i];
    atom->bond_atom[i][m] = partner;
    atom->bond_type[i][m] = bond_type;
    // Clear stale history data before reusing this slot for a new bond.
    if (n_histories > 0)
      for (auto &hist : histories)
        dynamic_cast<FixBondHistory *>(hist)->delete_history(i, m);
    atom->num_bond[i]++;
  };

  auto remove_special_12 = [&](int i, tagint partner) {
    tagint *slist = atom->special[i];
    for (int m = 0; m < atom->nspecial[i][0]; m++) {
      if (slist[m] != partner) continue;
      for (int n = m; n < atom->nspecial[i][2] - 1; n++) slist[n] = slist[n + 1];
      atom->nspecial[i][0]--;
      atom->nspecial[i][1]--;
      atom->nspecial[i][2]--;
      return true;
    }
    return false;
  };

  auto add_special_12 = [&](int i, tagint partner) {
    tagint *slist = atom->special[i];
    int n1 = atom->nspecial[i][0];
    int n2 = atom->nspecial[i][1];
    int n3 = atom->nspecial[i][2];

    // Match fix bond/create: if the new bonded partner already appears as a
    // 1-3 or 1-4 special neighbor.
    for (int m = n1; m < n3; m++) {
      if (slist[m] != partner) continue;
      for (int n = m; n < n3 - 1; n++) slist[n] = slist[n + 1];
      n3--;
      if (m < n2) n2--;
      break;
    }

    if (n3 == atom->maxspecial)
      error->one(FLERR, "New bond from fix {} exceeds special list size limit", style);
    for (int m = n3; m > n1; m--) slist[m] = slist[m - 1];
    slist[n1] = partner;
    atom->nspecial[i][0] = n1 + 1;
    atom->nspecial[i][1] = n2 + 1;
    atom->nspecial[i][2] = n3 + 1;
  };

  if (pivot >= 0 && pivot < atom->nlocal) {
    if (old_on_pivot && new_on_pivot) {
      if (!replace_one(pivot, tag_leave, tag_attack))
        error->one(FLERR, "Fix {} could not find old bond on pivot atom", style);
    } else {
      if (old_on_pivot && !delete_one(pivot, tag_leave))
        error->one(FLERR, "Fix {} could not find old bond on pivot atom", style);
      if (new_on_pivot) add_one(pivot, tag_attack);
    }
    if (!remove_special_12(pivot, tag_leave))
      error->one(FLERR, "Fix {} could not find old 1-2 special neighbor on pivot atom",
                 style);
    add_special_12(pivot, tag_attack);
  }

  if (leave >= 0 && leave < atom->nlocal) {
    if (old_on_leave && !delete_one(leave, tag_pivot))
      error->one(FLERR, "Fix {} could not find old bond on leaving atom", style);
    if (!remove_special_12(leave, tag_pivot))
      error->one(FLERR, "Fix {} could not find old 1-2 special neighbor on leaving atom",
                 style);
    atom->type[leave] = free_type;
  }

  if (attack >= 0 && attack < atom->nlocal) {
    if (new_on_attack) add_one(attack, tag_pivot);
    add_special_12(attack, tag_pivot);
    atom->type[attack] = bound_type;
  }
}

/* ----------------------------------------------------------------------
   Forward all owned atom type/bond/special data to ghost copies.  This is
   used once at the start of each MC invocation because normal atom border
   communication for molecular styles does not populate ghost bond and special
   topology.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::sync_topology_all()
{
  comm->forward_comm(this, topology_record_size());
}

/* ----------------------------------------------------------------------
   Forward the accepted triplet's updated owned topology to any ghost copies.
   Only the pivot, leaving A, and attacking A can change in an accepted move, so
   this accepted-move sync avoids a full fixed-width forward_comm pass.  Each
   owner packs its authoritative record, non-owners contribute zeros, and MPI_SUM
   recovers the unique owned records for all ranks to unpack into local/ghost copies.
   The sum operation preserves possible negative bond_type values.
   This keeps later trials in the same Ntry sequence from seeing stale
   free/bound labels or stale dynamic B-A bond/special entries.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::sync_topology_for_move(tagint tag_pivot, tagint tag_leave,
                                                  tagint tag_attack)
{
  const int nper = topology_record_size();
  const int nall = 3 * nper;
  tagint tags[3] = {tag_pivot, tag_leave, tag_attack};

  topology_send.assign(nall, 0);
  topology_recv.assign(nall, 0);

  for (int i = 0; i < 3; i++)
    pack_owned_topology_record(tags[i], &topology_send[i * nper]);

  MPI_Allreduce(topology_send.data(), topology_recv.data(), nall, MPI_LMP_TAGINT, MPI_SUM,
                world);

  for (int i = 0; i < 3; i++)
    unpack_topology_record(tags[i], &topology_recv[i * nper]);
}

/* ----------------------------------------------------------------------
   Fixed number of slots in one topology record.  The full ghost refresh packs
   this many double slots per communicated atom through forward_comm(); the
   accepted-move sync packs the same integer fields into tagint buffers for
   the three changed atom tags.
------------------------------------------------------------------------- */

int FixBondExchangeAssoc::topology_record_size() const
{
  return 5 + 2 * atom->bond_per_atom + atom->maxspecial;
}

/* ----------------------------------------------------------------------
   Pack the owner record for one changed atom.  Only the owning rank contributes
   nonzero data; other ranks leave the zero-filled record unchanged so MPI_SUM can
   recover the unique owner state without losing negative topology types.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::pack_owned_topology_record(tagint tag, tagint *buf) const
{
  int i = atom->map(tag);
  if (i < 0 || i >= atom->nlocal) return;

  int m = 0;
  buf[m++] = atom->type[i];
  buf[m++] = atom->num_bond[i];
  for (int k = 0; k < atom->bond_per_atom; k++) {
    buf[m++] = (k < atom->num_bond[i]) ? atom->bond_type[i][k] : 0;
    buf[m++] = (k < atom->num_bond[i]) ? atom->bond_atom[i][k] : 0;
  }
  buf[m++] = atom->nspecial[i][0];
  buf[m++] = atom->nspecial[i][1];
  buf[m++] = atom->nspecial[i][2];
  for (int k = 0; k < atom->maxspecial; k++)
    buf[m++] = (k < atom->nspecial[i][2]) ? atom->special[i][k] : 0;
}

/* ----------------------------------------------------------------------
   Unpack a compact owner record into a local or ghost copy of the changed atom.
   All ranks receive all three owner records, but only ranks with the atom in
   their current map need to store it.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::unpack_topology_record(tagint tag, const tagint *buf)
{
  int i = atom->map(tag);
  if (i < 0 || i >= atom->nlocal + atom->nghost) return;

  int m = 0;
  int atype = static_cast<int>(buf[m++]);
  if (atype <= 0)
    error->one(FLERR, "Fix {} received invalid compact topology record", style);
  atom->type[i] = atype;

  int nbond = static_cast<int>(buf[m++]);
  if (nbond > atom->bond_per_atom)
    error->one(FLERR, "Fix {} received too many compact ghost bonds", style);
  atom->num_bond[i] = nbond;
  for (int k = 0; k < atom->bond_per_atom; k++) {
    atom->bond_type[i][k] = static_cast<int>(buf[m++]);
    atom->bond_atom[i][k] = buf[m++];
  }

  int ns0 = static_cast<int>(buf[m++]);
  int ns1 = static_cast<int>(buf[m++]);
  int ns2 = static_cast<int>(buf[m++]);
  if (ns2 > atom->maxspecial)
    error->one(FLERR, "Fix {} received too many compact ghost special neighbors", style);
  atom->nspecial[i][0] = ns0;
  atom->nspecial[i][1] = ns1;
  atom->nspecial[i][2] = ns2;
  for (int k = 0; k < atom->maxspecial; k++)
    atom->special[i][k] = buf[m++];
}

/* ----------------------------------------------------------------------
   Pack owned type, bond list, and special list for the full ghost-topology
   refresh.  The accepted-move sync uses the compact tagint sync above; this
   fixed-width double record is retained for sync_topology_all().
------------------------------------------------------------------------- */

int FixBondExchangeAssoc::pack_forward_comm(int n, int *list, double *buf,
                                            int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;

  for (int ii = 0; ii < n; ii++) {
    int i = list[ii];

    buf[m++] = ubuf(atom->type[i]).d;
    buf[m++] = ubuf(atom->num_bond[i]).d;
    for (int k = 0; k < atom->bond_per_atom; k++) {
      int btype = (k < atom->num_bond[i]) ? atom->bond_type[i][k] : 0;
      tagint batom = (k < atom->num_bond[i]) ? atom->bond_atom[i][k] : 0;
      buf[m++] = ubuf(btype).d;
      buf[m++] = ubuf(batom).d;
    }

    buf[m++] = ubuf(atom->nspecial[i][0]).d;
    buf[m++] = ubuf(atom->nspecial[i][1]).d;
    buf[m++] = ubuf(atom->nspecial[i][2]).d;
    for (int k = 0; k < atom->maxspecial; k++) {
      tagint stag = (k < atom->nspecial[i][2]) ? atom->special[i][k] : 0;
      buf[m++] = ubuf(stag).d;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   Unpack owner topology into ghost atoms for the full ghost-topology refresh.
   The record layout matches pack_forward_comm() and contains no active flag.
------------------------------------------------------------------------- */

void FixBondExchangeAssoc::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;

  for (int i = first; i < last; i++) {
    atom->type[i] = (int) ubuf(buf[m++]).i;
    int nbond = (int) ubuf(buf[m++]).i;
    if (nbond > atom->bond_per_atom)
      error->one(FLERR, "Fix {} received too many ghost bonds", style);
    atom->num_bond[i] = nbond;
    for (int k = 0; k < atom->bond_per_atom; k++) {
      atom->bond_type[i][k] = (int) ubuf(buf[m++]).i;
      atom->bond_atom[i][k] = (tagint) ubuf(buf[m++]).i;
    }

    int ns0 = (int) ubuf(buf[m++]).i;
    int ns1 = (int) ubuf(buf[m++]).i;
    int ns2 = (int) ubuf(buf[m++]).i;
    if (ns2 > atom->maxspecial)
      error->one(FLERR, "Fix {} received too many ghost special neighbors", style);
    atom->nspecial[i][0] = ns0;
    atom->nspecial[i][1] = ns1;
    atom->nspecial[i][2] = ns2;
    for (int k = 0; k < atom->maxspecial; k++)
      atom->special[i][k] = (tagint) ubuf(buf[m++]).i;
  }
}

/* ----------------------------------------------------------------------
   Let users replace the internally created temperature compute, following
   the same "fix_modify ID temp compute-ID" convention as fix bond/swap.
------------------------------------------------------------------------- */

int FixBondExchangeAssoc::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "temp") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal fix_modify command");

    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete[] id_temp;

    id_temp = utils::strdup(arg[1]);
    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature)
      error->all(FLERR, "Could not find fix_modify temperature ID");
    if (temperature->tempflag == 0)
      error->all(FLERR, "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR, "Group for fix_modify temp != fix group");

    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   Global scalar output: accepted / accepted-plus-rejected Metropolis decisions.
   Trials skipped by pkin or with no proposal are tracked in the vector but are
   not Metropolis rejections.
------------------------------------------------------------------------- */

double FixBondExchangeAssoc::compute_scalar()
{
  double one[2], all[2];
  one[0] = static_cast<double>(naccept);
  one[1] = static_cast<double>(naccept + nreject);
  MPI_Allreduce(one, all, 2, MPI_DOUBLE, MPI_SUM, world);
  if (all[1] == 0.0) return 0.0;
  return all[0] / all[1];
}

/* ----------------------------------------------------------------------
   Global vector output for cumulative attempt, accept, reject, no-candidate,
   and pkin-skip counters.  The final vector entry returns the same acceptance
   ratio as compute_scalar().
------------------------------------------------------------------------- */

double FixBondExchangeAssoc::compute_vector(int n)
{
  double one, all;
  if (n == 0) one = static_cast<double>(ntrials);
  else if (n == 1) one = static_cast<double>(naccept);
  else if (n == 2) one = static_cast<double>(nreject);
  else if (n == 3) one = static_cast<double>(n_no_pivot);
  else if (n == 4) one = static_cast<double>(n_no_leaving);
  else if (n == 5) one = static_cast<double>(n_no_attacker);
  else if (n == 6) one = static_cast<double>(n_skip_pkin);
  else return compute_scalar();

  MPI_Allreduce(&one, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/* ----------------------------------------------------------------------
   Memory used by growable C++ work vectors.
------------------------------------------------------------------------- */

double FixBondExchangeAssoc::memory_usage()
{
  double bytes = 0.0;
  bytes += static_cast<double>(pivots.capacity()) * sizeof(int);
  bytes += static_cast<double>(attackers.capacity()) * sizeof(int);
  bytes += static_cast<double>(active_ranks.capacity()) * sizeof(int);
  bytes += static_cast<double>(topology_send.capacity()) * sizeof(tagint);
  bytes += static_cast<double>(topology_recv.capacity()) * sizeof(tagint);
  return bytes;
}
