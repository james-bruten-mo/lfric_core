##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

ROOT     = ../..
BIN_DIR  = $(ROOT)/bin

include $(ROOT)/make/include.mk

# Identify any scripts which need to be copied
SCRIPTS := $(shell find . ! -name '*~' ! -name Makefile )
SCRIPTS := $(patsubst ./%, $(BIN_DIR)/%, $(SCRIPTS) )

.PHONY: all
all: $(SCRIPTS)

# Copy scripts
$(BIN_DIR)/%: % | $(BIN_DIR)
	@echo -e $(VT_BOLD)Copying$(VT_RESET) $<
	cp $< $@

# Create directories

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
