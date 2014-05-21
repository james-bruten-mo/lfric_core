##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

.PHONY: test
test: all
	$(MAKE) -C src/test

.PHONY: all
all:
	$(MAKE) -C src/dynamo

.PHONY: run
run: test
	$(MAKE) -C run

.PHONY: doc docs
doc docs:
	$(MAKE) -C Docs

.PHONY: clean
clean:
	$(MAKE) -C src/dynamo clean
	$(MAKE) -C src/test clean
