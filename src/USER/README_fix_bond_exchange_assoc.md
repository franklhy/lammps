# fix bond/exchange/assoc

Prototype Markdown documentation for the USER package implementation in
`fix_bond_exchange_assoc.cpp`. The section order follows the LAMMPS manual page
style for fix commands.

## Syntax

```lammps
fix ID group-ID bond/exchange/assoc Nevery Ntry seed Rmin Rmax keyword values ...
```

- `ID`, `group-ID` are documented in the LAMMPS `fix` command.
- `bond/exchange/assoc` is the style name of this fix command.
- `Nevery` = attempt associative bond exchange every this many steps.
- `Ntry` = number of Monte Carlo trials each time the fix is invoked.
- `seed` = random number seed, a positive integer.
- `Rmin` = minimum old and new B-A distance for an exchange.
- `Rmax` = maximum old and new B-A distance for an exchange.
- Required keyword/value pairs: `pivot_type`, `free_type`, `bound_type`, `bond_type`, `pkin`.
- Optional keyword/value pair: `defensive_checks`.

Keyword syntax:

```text
pivot_type BTYPE
  BTYPE = atom type of the pivot B atom
free_type AFREE
  AFREE = atom type of the free A state
bound_type ABOUND
  ABOUND = atom type of the bound A state
bond_type BOND_TYPE
  BOND_TYPE = dynamic B-A bond type to exchange
pkin PKIN
  PKIN = kinetic prefactor, 0 < PKIN <= 1
defensive_checks yes|no
  enable additional consistency checks before committing accepted moves
```

## Examples

```lammps
fix xchg all bond/exchange/assoc 100 1000 12345 0.8 1.2 &
    pivot_type 3 free_type 4 bound_type 5 bond_type 2 pkin 1.0
```

This example attempts associative bond exchanges every 100 timesteps. Each
invocation performs up to 1000 MC trials. Atom type 3 is the pivot atom `B`, atom
type 4 is the free `A` state, atom type 5 is the bound `A` state, and bond type 2
is the dynamic `B-A` bond that is deleted and recreated.

## Description

This fix performs Monte Carlo associative bond exchange moves during a molecular
dynamics simulation. It is intended for dynamic networks where a pivot atom `B`
is initially bonded to A atom `Ai`, another A atom `Aj` is free, and the dynamic
bond can be transferred without changing the total number of dynamic bonds or
the total number of free and bound A states:

```text
Ai-B + Aj  ->  Ai + B-Aj
```

Here `Ai` is the leaving A atom: before the move it is bonded to `B` and has type
`bound_type`. `Aj` is the attacking A atom: before the move it is free and has
type `free_type`. In one accepted move, the old dynamic bond between `Ai` and
`B` is deleted, a new dynamic bond between `B` and `Aj` is created, and the two A
atoms exchange their state labels:

```text
delete bond Ai-B of type BOND_TYPE
create bond B-Aj of type BOND_TYPE
type[Ai] = AFREE
type[Aj] = ABOUND
```

Here `BTYPE` is the pivot atom type, `AFREE` is the free A atom type, `ABOUND` is
the bound A atom type, and `BOND_TYPE` is the dynamic B-A bond type. The three
atom type labels must be distinct. The `free_type` and `bound_type` atoms
represent the same chemical species in different bonding states; they are
required to have identical masses and identical pair interactions with every
atom type.

A check for possible exchanges is performed every `Nevery` timesteps. Each
invocation performs `Ntry` trials. A trial first applies the move-independent
kinetic gate `pkin`. If a uniform random number is larger than `pkin`, the trial
is skipped before selecting a pivot or building candidates. Otherwise the fix
selects one pivot atom uniformly from all owned atoms in the fix group with type
`pivot_type`. The active MPI rank is the rank that owns that pivot.

For the selected pivot `B`, the leaving atom `Ai` must satisfy all of the
following criteria:

- `Ai` is bonded to `B` by a bond of type `bond_type`.
- `Ai` has type `bound_type`.
- `Ai` is in the fix group.
- The old distance `r_BAi` is in the interval `[Rmin,Rmax]`.

The attacking atom `Aj` must satisfy all of the following criteria:

- `Aj` has type `free_type`.
- `Aj` is in the fix group.
- `Aj` is in the neighbor list of `B`.
- The new distance `r_BAj` is in the interval `[Rmin,Rmax]`.

The same symmetric reactive window is applied to the old and new B-A distances.
This is required so the reverse move is proposal-accessible immediately after an
accepted forward move. Do not choose `Rmax` so close to a bond singularity, such
as the FENE maximum extension, that trial energy evaluations themselves become
problematic.

If there are multiple attacking atoms, one is chosen uniformly. In the current
implementation, both A and B are assumed to be monovalent with respect to the
dynamic `bond_type`. Thus a `free_type` atom should have no direct B-A bond, a
`bound_type` atom should have exactly one dynamic B-A bond, and a `pivot_type`
atom should have at most one dynamic B-A bond. The fix performs a one-time
startup topology check for these invariants before the first MC attempt.

For a generated proposal, the energy difference `dU` is computed locally on the
active rank. Here `dU` means the potential energy of the proposed state minus
the potential energy of the current state:

```text
dU = U_after - U_before
```

For the current symmetric implementation, this is evaluated as:

```text
dU = [U_bond(B,Aj) - U_bond(B,Ai)]
   + [U_pair^0(B,Ai) + U_pair^1(B,Aj)
      - U_pair^1(B,Ai) - U_pair^0(B,Aj)]
```

`U_pair^0` is the normal, unbonded pair energy, and `U_pair^1` is the pair energy
evaluated with the 1-2 `special_bonds` scaling for a directly bonded pair. Pair
terms outside the pair cutoff contribute zero. In this symmetric implementation,
ordinary pair interactions involving `free_type` and `bound_type` A atoms cancel
because those two atom types are required to have identical pair interactions
with every atom type. For `pair_style lj/cut`, the fix checks equality of the
extracted epsilon, sigma, and cutoff values. For other supported two-body pair
styles, the fix compares sampled calls to `Pair::single()` at `Rmin`,
`(Rmin+Rmax)/2`, and `Rmax`.

The acceptance probability for trials that pass the `pkin` gate and generate a
proposal is:

```text
Pacc = min(1, exp(-dU/(kB*T)))
```

where `kB` is the Boltzmann constant for the current LAMMPS unit style and `T` is
the current temperature computed by this fix. Since `pkin` is applied before
proposal generation and does not depend on the proposed move, the total
acceptance probability is equivalent to `pkin` times the Metropolis factor above.
In a more general Metropolis-Hastings move, the acceptance probability would also
include a proposal-probability ratio `qrev/qfwd`, where `qfwd` is the probability
of proposing the current forward move and `qrev` is the probability of proposing
the exact reverse move from the proposed state. For this simple monovalent
symmetric exchange, the forward and reverse proposal counts are the same, so
`qrev/qfwd = 1` and the ratio does not appear explicitly in `Pacc`.

This fix computes a temperature each time it is invoked for use by the Boltzmann
criterion. To do this, the fix creates its own compute of style `temp`, as if
this command had been issued:

```lammps
compute fix-ID_temp all temp
```

The ID of the new compute is the fix ID with underscore + `temp` appended. You
can use `fix_modify` with the `temp` keyword to assign another temperature
compute to this fix:

```lammps
fix_modify xchg temp myTemp
```

When an accepted move changes the topology, the fix updates atom types, bond
lists, and the 1-2 special-neighbor lists for the affected atoms. The accepted
atom tags are broadcast to all MPI ranks, each rank modifies only atoms it owns,
and the updated owned topology is then communicated to ghost copies. Thus no
reverse communication of topology changes from ghosts to owners is needed.

The initial Bonds section does not need to list the pivot atom first for dynamic
B-A bonds. The fix finds existing dynamic bonds from the 1-2 special-neighbor
list and removes the old bond from whichever endpoint actually stores it. New
dynamic bonds created by this fix follow the same convention as `fix bond/create`:
with `newton bond on`, the lower-tag endpoint stores the new bond;
with `newton bond off`, both endpoints store it.

The optional `defensive_checks` keyword enables an additional consistency check
before committing an accepted move. With this option enabled, the fix verifies
that the old dynamic bond has the expected number of stored bond-list entries
before deleting it. This option is useful while validating new systems or
debugging topology issues, but it adds extra communication and is off by default.

You can dump snapshots of the current bond topology via `dump local`, for
example with `compute property/local btype batom1 batom2`.

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to binary restart files. Topology and
atom type changes that have already been accepted are stored as ordinary atom
topology in data and restart files, but the random number state and cumulative
fix statistics are not stored. Thus restarted simulations are not exact bitwise
continuations of the Monte Carlo sequence, although they should produce
equivalent statistical behavior.

The `fix_modify temp` option is supported by this fix. You can use it to assign a
temperature compute for the Boltzmann criterion.

This fix computes a global scalar and a global vector of length 8, which can be
accessed by LAMMPS output commands. The scalar and vector values calculated by
this fix are intensive.

The global scalar is the cumulative Metropolis acceptance ratio for trials that
passed the `pkin` gate and generated a proposal. Trials skipped by the `pkin`
gate or trials without a valid proposal are reported separately in the vector.

The vector values are:

1. `f_ID[1]`: cumulative MC trials
2. `f_ID[2]`: cumulative accepted moves
3. `f_ID[3]`: cumulative rejected proposed moves
4. `f_ID[4]`: cumulative no-pivot trials
5. `f_ID[5]`: cumulative no-leaving-atom trials
6. `f_ID[6]`: cumulative no-attacking-atom trials
7. `f_ID[7]`: cumulative trials skipped by the `pkin` gate
8. `f_ID[8]`: same Metropolis acceptance ratio as the global scalar

No parameter of this fix can be used with the `start/stop` keywords of the `run`
command. This fix is not invoked during energy minimization.

## Restrictions

This fix is part of the USER package. It is only enabled if LAMMPS is built with
that package.

This fix requires molecular atom styles with bonds, atom IDs, an atom map, a pair
style, and a bond style.

This implementation supports only the monovalent symmetric A-state-transfer
model described above. It does not support systems with angles, dihedrals, or
impropers involving the dynamic bond. It does not support kspace/long-range
electrostatics or many-body pair styles. It requires pair styles that provide a
usable `Pair::single()` energy calculation.

The `free_type` and `bound_type` atom types must have identical masses and
identical pair interactions with every atom type. The 1-3 and 1-4 weights of the
`special_bonds` command must be 1.0.

## Related commands

`fix bond/create`, `fix bond/break`, `fix bond/react`, `fix bond/swap`, `dump
local`, `compute property/local`, `special_bonds`, `fix_modify`

## Default

The option default is:

```text
defensive_checks = no
```

## Validation Examples

Example input/data files and the completed validation checklist are in:

```text
examples/PACKAGES/user/README_fix_bond_exchange_assoc_tests.md
```
