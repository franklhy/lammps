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
  YAFF)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(BUILD_TOOLS ON CACHE BOOL "" FORCE)

### GPU package, from gpu-cuda.cmake
set(PKG_GPU ON CACHE BOOL "Build GPU package" FORCE)
set(GPU_API "cuda" CACHE STRING "APU used by GPU package" FORCE)
set(GPU_PREC "mixed" CACHE STRING "" FORCE)

set(CUDA_NVCC_FLAGS "-allow-unsupported-compiler" CACHE STRING "" FORCE)
set(CUDA_NVCC_FLAGS_DEBUG "-allow-unsupported-compiler" CACHE STRING "" FORCE)
set(CUDA_NVCC_FLAGS_MINSIZEREL "-allow-unsupported-compiler" CACHE STRING "" FORCE)
set(CUDA_NVCC_FLAGS_RELWITHDEBINFO "-allow-unsupported-compiler" CACHE STRING "" FORCE)
set(CUDA_NVCC_FLAGS_RELEASE "-allow-unsupported-compiler" CACHE STRING "" FORCE)


### KOKKS package, from kokkos-cuda.cmake
set(PKG_KOKKOS ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA   ON CACHE BOOL "" FORCE)
set(BUILD_OMP ON CACHE BOOL "" FORCE)
get_filename_component(NVCC_WRAPPER_CMD ${CMAKE_CURRENT_SOURCE_DIR}/../lib/kokkos/bin/nvcc_wrapper ABSOLUTE)
set(CMAKE_CXX_COMPILER ${NVCC_WRAPPER_CMD} CACHE FILEPATH "" FORCE)

# If KSPACE is also enabled, use CUFFT for FFTs
set(FFT_KOKKOS "CUFFT" CACHE STRING "" FORCE)

# hide deprecation warnings temporarily for stable release
set(Kokkos_ENABLE_DEPRECATION_WARNINGS OFF CACHE BOOL "" FORCE)
