##############################################################################
# (c) Crown copyright 2018 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# This file is for any project specific build settings to be applied
# via the Makefile.

$(info Coupled miniapp project specials)

export PRE_PROCESS_MACROS += MCT
export EXTERNAL_STATIC_LIBRARIES += psmile.MPI1 mct mpeu scrip
