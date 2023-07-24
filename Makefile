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

# --------------------------------------------------------
# Program build logic

# Programs to build, assuming each has a corresponding *.chpl file
PROGRAMS = \
  diffusion/heat_2D \
  diffusion/heat_2D_dist \
  diffusion/heat_2D_dist_stencil \
  diffusion/heat_2D_dist_buffers \
  diffusion/heat_2D_dist_exchanges \
  diffusion/heat_2D_dist_exchanges_abstracted \
  image_analysis/main \
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

# The rule for building any example.
%: %.chpl
	$(CHPL) $@.chpl --fast

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
	./$<
	-@if test -f $(OUTPUT) ; then \
	     echo $(OUTPUT) : ; cat $(OUTPUT) ; \
	  fi

clean:
	rm -f $(PROGRAMS) $(OUTPUT)

force:

.PHONY: clean all run force
