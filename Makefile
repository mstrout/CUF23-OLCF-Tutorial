# This Makefile demonstrates the recommended way to build small Chapel programs.
# Note this uses some GNU make extensions for conciseness.
#
# To run an example like heat_2D do the following command:
#    make run-heat_2D
# If you just want to compile then do
#    make heat_2D
#

CHPL=chpl

ifneq ($(strip $(SLURM_JOB_NUM_NODES)),)
  # Inside an salloc session
  EXECARG ?= -nl $(SLURM_JOB_NUM_NODES)
else ifneq (, $(shell which srun 2> /dev/null))
  # From a login node with SLURM
  EXECARG ?= -nl 1
else
  # Unknown system (e.g. AWS or other)
  EXECARG ?= -nl 1
endif


# --------------------------------------------------------
# Program build logic

# Programs to build, assuming each has a corresponding *.chpl file
PROGRAMS = \
  basics-coforall \
  basics-distarr \
  basics-for \
  basics-on \
  heat_2D \
  heat_2D_dist \
  heat_2D_dist_stencil \
  heat_2D_dist_buffers \
  heat_2D_dist_exchanges \
  heat_2D_dist_exchanges_abstracted \
  heat_1D \
  heat_1D_dist \
  heat_1D_tasks \
  gpuExample \
  hello-dist-node-names \
  hello \
  hello6-taskpar-dist \
  hellopar \
  kmer \
  parfilekmer \
  stream-ep \
  writelnExamples

all: $(PROGRAMS)

EXTRA_FLAGS ?= 

# The rule for building any example.
%: %.chpl
	$(CHPL) $@.chpl -o $@ --fast --no-warnings $(EXTRA_FLAGS)

# --------------------------------------------------------
# Everything below is convenience targets for usability

# runs all the programs currently built in the working directory
run:
	@for f in $(PROGRAMS) ; do \
	  if test -x $$f ; then \
	    ( set -x ; ./$$f $(EXECARG) ; ) ; \
	  fi ; \
	done

# builds and runs a particular test (eg `make run-ex0`)
run-%: % force
	./$< $(EXECARG)

clean:
	rm -f $(PROGRAMS) $(PROGRAMS:=_real) main main_real

force:

.PHONY: clean all run force
