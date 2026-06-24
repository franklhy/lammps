.. index:: fix bond/exchange/assoc

fix bond/exchange/assoc command
================================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID bond/exchange/assoc Nevery Ntry seed Rmin Rmax keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* bond/exchange/assoc = style name of this fix command
* Nevery = attempt associative bond exchange every this many steps
* Ntry = number of Monte Carlo trials each time the fix is invoked
* seed = random number seed (positive integer)
* Rmin = minimum old and new B-A distance for an exchange (distance units)
* Rmax = maximum old and new B-A distance for an exchange (distance units)
* required keyword/value pairs = *pivot_type*, *free_type*, *bound_type*, *bond_type*, *pkin*
* optional keyword/value pair = *defensive_checks*

  .. parsed-literal::

       *pivot_type* value = BTYPE
         BTYPE = atom type of the pivot B atom (integer or type label)
       *free_type* value = AFREE
         AFREE = atom type of the free A state (integer or type label)
       *bound_type* value = ABOUND
         ABOUND = atom type of the bound A state (integer or type label)
       *bond_type* value = BOND_TYPE
         BOND_TYPE = dynamic B-A bond type to exchange (integer or type label)
       *pkin* value = PKIN
         PKIN = kinetic prefactor, 0 < PKIN <= 1
       *defensive_checks* value = *yes* or *no*
         enable additional consistency checks before committing accepted moves

Examples
""""""""

.. code-block:: LAMMPS

   fix xchg all bond/exchange/assoc 100 1000 12345 0.8 1.2 &
       pivot_type 3 free_type 4 bound_type 5 bond_type 2 pkin 1.0

Description
"""""""""""

Perform Monte Carlo associative bond exchange moves during a molecular dynamics
simulation.  This fix is intended for dynamic networks where a pivot atom
:math:`B` is initially bonded to A atom ``Ai``, another A atom ``Aj`` is free,
and the dynamic bond can be transferred without changing the total number of
dynamic bonds or the total number of free and bound A states:

.. parsed-literal::

   Ai-B + Aj  ->  Ai + B-Aj

Here ``Ai`` is the leaving A atom: before the move it is bonded to :math:`B` and
has type *bound_type*.  ``Aj`` is the attacking A atom: before the move it is
free and has type *free_type*.  In one accepted move, the old dynamic bond
between ``Ai`` and :math:`B` is deleted, a new dynamic bond between :math:`B` and
``Aj`` is created, and the two A atoms exchange their state labels:

.. parsed-literal::

   delete bond Ai-B of type BOND_TYPE
   create bond B-Aj of type BOND_TYPE
   type[Ai] = AFREE
   type[Aj] = ABOUND

Here *BTYPE* is the pivot atom type, *AFREE* is the free A atom type, *ABOUND*
is the bound A atom type, and *BOND_TYPE* is the dynamic B-A bond type.  The
three atom type labels must be distinct.  The *free_type* and *bound_type* atoms
represent the same chemical species in different bonding states; they are
required to have identical masses and identical pair interactions with every
atom type.

A check for possible exchanges is performed every *Nevery* timesteps.  Each
invocation performs *Ntry* trials.  A trial first applies the move-independent
kinetic gate *pkin*.  If a uniform random number is larger than *pkin*, the trial
is skipped before selecting a pivot or building candidates.  Otherwise the fix
selects one pivot atom uniformly from all owned atoms in the fix group with type
*pivot_type*.  The active MPI rank is the rank that owns that pivot.

For the selected pivot :math:`B`, the leaving atom ``Ai`` must satisfy all of
the following criteria:

* ``Ai`` is bonded to :math:`B` by a bond of type *bond_type*
* ``Ai`` has type *bound_type*
* ``Ai`` is in the fix group
* the old distance ``r_BAi`` is in the interval [*Rmin*, *Rmax*]

The attacking atom ``Aj`` must satisfy all of the following criteria:

* ``Aj`` has type *free_type*
* ``Aj`` is in the fix group
* ``Aj`` is in the neighbor list of :math:`B`
* the new distance ``r_BAj`` is in the interval [*Rmin*, *Rmax*]

The same symmetric reactive window is applied to the old and new B-A distances.
This is required so the reverse move is proposal-accessible immediately after an
accepted forward move.  Do not choose *Rmax* so close to a bond singularity, such
as the FENE maximum extension, that trial energy evaluations themselves become
problematic.

If there are multiple attacking atoms, one is chosen uniformly.  In the current
implementation, both A and B are assumed to be monovalent with respect to the
dynamic *bond_type*.  Thus a *free_type* atom should have no direct B-A bond, a
*bound_type* atom should have exactly one dynamic B-A bond, and a *pivot_type*
atom should have at most one dynamic B-A bond.  The fix performs a one-time
startup topology check for these invariants before the first MC attempt.

For a generated proposal, the energy difference :math:`\Delta U` is computed
locally on the active rank.  The calculation includes the exchanged bond energy
and the pair-energy change associated with the old ``Ai-B`` and new ``B-Aj``
direct pairs.  In this symmetric implementation, ordinary pair interactions
involving *free_type* and *bound_type* A atoms cancel because those two atom
types are required to have identical pair interactions with every atom type.  For
*pair_style lj/cut*, the
fix checks equality of the extracted epsilon, sigma, and cutoff values.  For
other supported two-body pair styles, the fix compares sampled calls to
``Pair::single()`` at *Rmin*, :math:`(Rmin+Rmax)/2`, and *Rmax*.

The acceptance probability for trials that pass the *pkin* gate and generate a
proposal is

.. math::

   P_\mathrm{acc} = \min\left[1, \exp\left(-\Delta U / k_B T\right)\right]

where :math:`T` is the current temperature computed by this fix.  Since *pkin*
is applied before proposal generation and does not depend on the proposed move,
the total acceptance probability is equivalent to *pkin* times the Metropolis
factor above.  The simple monovalent symmetric proposal has
:math:`q_\mathrm{rev}/q_\mathrm{fwd} = 1`.

This fix computes a temperature each time it is invoked for use by the
Boltzmann criterion.  To do this, the fix creates its own compute of style
*temp*, as if this command had been issued:

.. code-block:: LAMMPS

   compute fix-ID_temp all temp

The ID of the new compute is the fix ID with underscore + "temp" appended.  You
can use :doc:`fix_modify <fix_modify>` with the *temp* keyword to assign another
temperature compute to this fix.

When an accepted move changes the topology, the fix updates atom types, bond
lists, and the 1-2 special-neighbor lists for the affected atoms.  The accepted
atom tags are broadcast to all MPI ranks, each rank modifies only atoms it owns,
and the updated owned topology is then communicated to ghost copies.  Thus no
reverse communication of topology changes from ghosts to owners is needed.

The initial Bonds section does not need to list the pivot atom first for dynamic
B-A bonds.  The fix finds existing dynamic bonds from the 1-2 special-neighbor
list and removes the old bond from whichever endpoint actually stores it.  New
dynamic bonds created by this fix follow the same convention as
:doc:`fix bond/create <fix_bond_create>`: with ``newton bond on``, the lower-tag
endpoint stores the new bond; with ``newton bond off``, both endpoints store it.

The optional *defensive_checks* keyword enables an additional consistency check
before committing an accepted move.  With this option enabled, the fix verifies
that the old dynamic bond has the expected number of stored bond-list entries
before deleting it.  This option is useful while validating new systems or
debugging topology issues, but it adds extra communication and is off by default.

You can dump snapshots of the current bond topology via the :doc:`dump local <dump>`
command, for example with :doc:`compute property/local <compute_property_local>`
using the *btype*, *batom1*, and *batom2* attributes.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
Topology and atom type changes that have already been accepted are stored as
ordinary atom topology in data and restart files, but the random number state and
cumulative fix statistics are not stored.  Thus restarted simulations are not
exact bitwise continuations of the Monte Carlo sequence, although they should
produce equivalent statistical behavior.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this fix.  You
can use it to assign a temperature compute for the Boltzmann criterion.

This fix computes a global scalar and a global vector of length 8, which can be
accessed by various :doc:`output commands <Howto_output>`.  The scalar and vector
values calculated by this fix are "intensive".

The global scalar is the cumulative Metropolis acceptance ratio for trials that
passed the *pkin* gate and generated a proposal.  Trials skipped by the *pkin*
gate or trials without a valid proposal are reported separately in the vector.

The vector values are:

  (1) cumulative MC trials
  (2) cumulative accepted moves
  (3) cumulative rejected proposed moves
  (4) cumulative no-pivot trials
  (5) cumulative no-leaving-atom trials
  (6) cumulative no-attacking-atom trials
  (7) cumulative trials skipped by the *pkin* gate
  (8) same Metropolis acceptance ratio as the global scalar

No parameter of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the USER package.  It is only enabled if LAMMPS was built with
that package.  See the :doc:`Build package <Build_package>` doc page for more
info.

This fix requires molecular atom styles with bonds, atom IDs, an atom map, a
pair style, and a bond style.

This implementation supports only the monovalent symmetric A-state-transfer
model described above.  It does not support systems with angles, dihedrals, or
impropers involving the dynamic bond.  It does not support kspace/long-range
electrostatics or many-body pair styles.  It requires pair styles that provide a
usable ``Pair::single()`` energy calculation.

The *free_type* and *bound_type* atom types must have identical masses and
identical pair interactions with every atom type.  The 1-3 and 1-4 weights of
the :doc:`special_bonds <special_bonds>` command must be 1.0.

Related commands
""""""""""""""""

:doc:`fix bond/create <fix_bond_create>`,
:doc:`fix bond/break <fix_bond_break>`,
:doc:`fix bond/react <fix_bond_react>`,
:doc:`fix bond/swap <fix_bond_swap>`, :doc:`dump local <dump>`,
:doc:`compute property/local <compute_property_local>`,
:doc:`special_bonds <special_bonds>`, :doc:`fix_modify <fix_modify>`

Default
"""""""

The option default is defensive_checks = no.
