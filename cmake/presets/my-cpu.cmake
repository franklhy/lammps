# preset that turns on a wide range of packages, some of which require
# external libraries. Compared to all_on.cmake some more unusual packages
# are removed. The resulting binary should be able to run most inputs.

set(ALL_PACKAGES
  ASPHERE
  BODY
  BROWNIAN
  CLASS2
  COLLOID
  DIELECTRIC
  DIPOLE
  EXTRA-COMMAND
  EXTRA-COMPUTE
  EXTRA-DUMP
  EXTRA-FIX
  EXTRA-MOLECULE
  EXTRA-PAIR
  KSPACE
  MANYBODY
  MC
  MISC
  MOLECULE
  OPENMP
  OPT
  ORIENT
  REACTION
  REAXFF
  REPLICA
  RIGID
  USER)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(BUILD_TOOLS ON CACHE BOOL "" FORCE)
