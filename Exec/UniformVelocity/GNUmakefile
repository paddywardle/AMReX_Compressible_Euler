# This should point to the location of the amrex source code.
AMREX_HOME = ../../../amrex_20_11/amrex

USE_EB =FALSE
PRECISION  = DOUBLE
PROFILE    = FALSE

# Debug mode compiles without optimisation, so a debugger returns
# accurate line information, without optimising out variables.
# Additionally, additional bounds checking is implemented in AMReX's
# debug mode
DEBUG      = TRUE
#DEBUG      = FALSE

# Dimensionality is specified at compile time - this example has been written for 2D only
DIM        = 2

# Compiler can be specified here
COMP	   = gnu

# MPI and OpenMP compatibility must be included at compile time.  OpenMP will also require appropriate pragmas.
USE_MPI    = TRUE
USE_OMP    = FALSE

Bpack   := ./Make.package 
Blocs   := . 

include ../Make.Adv
