#
# First build all of the CUDA projects, producing a series of algorithm
# libraries. Then build and link the main project against these.
#

ALGORITHMS := $(shell find src/cuda -name Makefile)
PROJECTS   := $(ALGORITHMS)
CABAL      := cabal

ifeq ($(dbg),1)
    CABALFLAGSCONF += -fdebug
endif
ifeq ($(verbose),1)
    CABALFLAGSBUILD += --verbose
    VERBOSE    :=
else
    VERBOSE    := @
endif


%.do :
	$(MAKE) -C $(dir $*) $(MAKECMDGOALS)

all : $(addsuffix .do,$(PROJECTS))
	$(VERBOSE)$(CABAL) configure $(CABALFLAGSCONF)
	$(VERBOSE)$(CABAL) build $(CABALFLAGSBUILD)
	@echo "Finished building all"

clean : $(addsuffix .do,$(PROJECTS))
	$(VERBOSE)$(CABAL) clean
	@echo "Finished cleaning all"

clobber : $(addsuffix .do,$(PROJECTS))
	$(VERBOSE)$(CABAL) clean
	@echo "Finished cleaning all"

cudabrush:
	rm ./dist/build/hfx/hfx-tmp/Foreign/CUDA/Algorithms.o
	rm ./dist/build/hfx/hfx-tmp/Foreign/CUDA/Util.o
