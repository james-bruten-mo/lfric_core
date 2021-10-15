##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the Intel Fortran compiler.
##############################################################################

$(info Project specials for Intel compiler)

export FFLAGS_UM_PHYSICS = -r8
# Remove -check-all option as it causes very slow runs due to a lot of array
# temporary warnings caused by UM code 
FFLAGS_RUNTIME            = -fpe0
FFLAGS_FORTRAN_STANDARD   = 

# NOTE: The -qoverride-limits option contained in $(FFLAGS_INTEL_FIX_ARG) is
# not currently applied here. This is a temporary workaround for #2340 in
# which it was found to be inadvertently applied to UM code preventing
# compilation of a UKCA module
#$(info LFRic compile options required for files with OpenMP - see Ticket 1490)
#%psy.o %psy.mod:   export FFLAGS += $(FFLAGS_INTEL_FIX_ARG)
#psy/%.o psy/%.mod: export FFLAGS += $(FFLAGS_INTEL_FIX_ARG)
