# This Makefile demonstrates the recommended way to build simple Chapel programs.
# Note this uses some GNU make extensions for conciseness.
#
# To use this makefile, set the CHPL_INSTALL variable to the chpl install directory, e.g.
# make CHPL_INSTALL=<myinstalldir> hello-world
# or (for bash)
# export CHPL_INSTALL=<myinstalldir>; make hello-world

ifeq ($(CHPL_INSTALL),)
$(warning CHPL_INSTALL environment variable is not set, assuming chpl is in the PATH)
CHPL=chpl
else
CHPL=$(CHPL_INSTALL)/bin/chpl
endif

ifeq ($(CHPL_COMM), "none")
EXECARG = ""
else
EXECARG = "-nl1"
endif

# --------------------------------------------------------
# Program build logic

# Programs to build, assuming each has a corresponding *.chpl file
PROGRAMS = \
  heat_2D \
  heat_2D_dist \
  heat_2D_dist_stencil \
  heat_2D_dist_buffers \
  heat_2D_dist_exchanges \
  heat_2D_dist_exchanges_abstracted \
  heat_1D \
  heat_1D_dist \
  heat_1D_tasks \
  main \
  gpuExample \
  hello-dist-node-names \
  hello \
  hello6-taskpar-dist \
  hellopar \
  kmer \
  parfilekmer \
  stream-ep \
  writelnExamples

REALS = \
  heat_2D_real \
  heat_2D_dist_real \
  heat_2D_dist_stencil_real \
  heat_2D_dist_buffers_real \
  heat_2D_dist_exchanges_real \
  heat_2D_dist_exchanges_abstracted_real \
  heat_1D_real \
  heat_1D_dist_real \
  heat_1D_tasks_real \
  main_real \
  gpuExample_real \
  hello-dist-node-names_real \
  hello_real \
  hello6-taskpar-dist_real \
  hellopar_real \
  kmer_real \
  parfilekmer_real \
  stream-ep_real \
  writelnExamples_real

all: $(PROGRAMS)

# The rule for building any example.
%: %.chpl
	$(CHPL) $@.chpl --fast --no-warnings

# --------------------------------------------------------
# Everything below is convenience targets for usability

OUTPUT = output.txt

# runs all the programs currently built in the working directory
run:
	@for f in $(PROGRAMS) ; do \
	  if test -x $$f ; then \
	    rm -f $(OUTPUT) ; \
	    ( set -x ; ./$$f ; ) ; \
	    if test -f $(OUTPUT) ; then \
	      echo $(OUTPUT) : ; cat $(OUTPUT) ; \
	    fi ; \
	  fi ; \
	done

# builds and runs a particular test (eg `make run-ex0`)
run-%: % force
	-@rm -f $(OUTPUT)
	./$< $(EXECARG)
	-@if test -f $(OUTPUT) ; then \
	     echo $(OUTPUT) : ; cat $(OUTPUT) ; \
	  fi

clean:
	rm -f $(PROGRAMS) $(OUTPUT) $(REALS)

force:

.PHONY: clean all run force
