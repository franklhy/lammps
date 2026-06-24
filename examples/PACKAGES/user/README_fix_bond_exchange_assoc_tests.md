# Tests for fix bond/exchange/assoc

This directory contains validation inputs for the associative bond-exchange fix.
All tests listed below have passed in local testing.

The tests are small, explicit LAMMPS input/data files rather than an automated
test harness.  Most inputs should be run from their own subdirectory so relative
data, dump, log, and restart filenames resolve cleanly.

## Validation Status

| Step | Test | Status |
| --- | --- | --- |
| 1 | Syntax/build check | passed |
| 2 | Minimal symmetric three-atom exchange | passed |
| 3 | Bond-energy acceptance | passed |
| 4 | `pkin` scaling | passed |
| 5 | Long-run conservation | passed |
| 6 | MPI cross-rank exchange sanity | passed |
| 7 | Defensive-check failures | passed |
| 8 | Data-file and binary-restart continuation | passed |

## 1. Syntax/build Check

Run a syntax-only compile check for `USER/fix_bond_exchange_assoc.cpp`, or the
normal CPU build when convenient.

Expected result: the fix compiles without errors.

## 2. Minimal Symmetric Three-Atom Exchange

Directory: `minimal_3_atom_exchange/`

Files:

- `in.fix_bond_exchange_assoc_minimal`
- `data.fix_bond_exchange_assoc_minimal`

System: one pivot `B`, one `A_bound`, and one `A_free`.  The old and new `B-A`
distances are symmetric, pair epsilon is zero, and `pkin = 1`.

Expected result:

- acceptance is near 1 for zero energy difference;
- exactly one dynamic `B-A` bond exists at all times;
- `A_free` and `A_bound` counts are conserved;
- no atom violates the monovalent dynamic-bond assumption.

## 3. Bond-Energy Acceptance

Directory: `bond_energy_acceptance/`

Files:

- `in.fix_bond_exchange_assoc_bond_energy`
- `data.fix_bond_exchange_assoc_bond_energy`

System: one `B`, one `A_bound`, and one `A_free`.  Pair interactions are zeroed
so the exchange probability is controlled by harmonic bond energy.  The uphill
proposal has `r_old = 1.0`, `r_new = 1.1`, harmonic `K = 30.0`, and `dU = 0.3`.

Expected result:

- uphill acceptance probability is `p = exp(-0.3/Tcurrent)`;
- the downhill reverse move is accepted with probability 1;
- for the fixed two-state geometry, the cumulative acceptance ratio `f_xchg`
  approaches `2*p/(1+p)`.

## 4. pkin Scaling

Directory: `pkin_scaling/`

Files:

- `in.fix_bond_exchange_assoc_pkin_1.0`
- `in.fix_bond_exchange_assoc_pkin_0.1`
- `in.fix_bond_exchange_assoc_pkin_0.01`
- `data.fix_bond_exchange_assoc_pkin`

System: the same symmetric zero-energy exchange is run with different `pkin`
values.

Expected result:

- Metropolis acceptance among generated proposals remains near 1 for zero `dU`;
- cumulative accepted moves `f_xchg[2]` scale approximately with `pkin`;
- cumulative skipped trials `f_xchg[7]` scale approximately with `1-pkin`;
- runtime decreases as `pkin` decreases.

Each input writes a distinct log file so the three runs can be launched from the
same directory without overwriting each other.

## 5. Long-Run Conservation

Directory: `long_run_conservation/`

Files:

- `in.fix_bond_exchange_assoc_conservation`
- `data.fix_bond_exchange_assoc_conservation`

System: four independent `B/A_bound/A_free` triplets.  Each triplet has two
symmetric possible `B-A` distances.

Expected conserved values:

- `c_nfree = 4`
- `c_nbound = 4`
- `c_ndynamic = 4`

The input uses `thermo_modify norm no` so these conserved quantities are printed
as counts rather than per-atom fractions.

## 6. MPI Cross-Rank Exchange Sanity

Directory: `mpi_cross_rank_exchange/`

Files:

- `in.fix_bond_exchange_assoc_mpi_cross_rank`
- `data.fix_bond_exchange_assoc_mpi_cross_rank`

Run with two MPI ranks.  The input uses `processors 2 1 1`, puts the pivot just
left of the x-domain split, and puts the attacking free A just right of the
split.

Expected result:

- no crashes or stale ghost-topology errors;
- the atom dump shows the pivot and attacking A owned by different ranks in the
  initial configuration;
- conserved counts match;
- acceptance statistics are plausible.

The atom dump includes `compute property/atom proc` as `c_owner`, so ownership is
visible directly in the dump.

## 7. Defensive-Check Failures

Directory: `defensive_checks/`

These are expected-failure tests.  A test passes when LAMMPS exits early with the
documented error substring.

Inputs:

- `in.fix_bond_exchange_assoc_missing_bound_bond`
  Expected error: `requires each bound_type atom to have one dynamic B-A bond`
- `in.fix_bond_exchange_assoc_free_direct_bond`
  Expected error: `found a free_type atom directly bonded to a pivot_type atom`
- `in.fix_bond_exchange_assoc_overvalent_pivot`
  Expected error: `requires each pivot_type atom to have at most one dynamic B-A bond`

## 8. Restart/Topology Robustness

Directory: `restart_topology/`

Files:

- `in.fix_bond_exchange_assoc_restart_topology`
- `in.fix_bond_exchange_assoc_restart_file`
- `data.fix_bond_exchange_assoc_restart_topology`

`in.fix_bond_exchange_assoc_restart_topology` is the data-file continuation
test.  It runs exchanges, writes
`data.fix_bond_exchange_assoc_restart_topology.stage1`, clears LAMMPS, reads the
generated data file, recreates the force field and fix, and continues exchanging.

`in.fix_bond_exchange_assoc_restart_file` is the binary-restart continuation
test.  It runs exchanges, writes
`restart.fix_bond_exchange_assoc_restart_file.stage1`, clears LAMMPS, reads the
restart file, recreates the computes and fix, and continues exchanging.

Expected result:

- topology remains valid after continuation;
- `A_free`, `A_bound`, and dynamic-bond counts remain conserved.
