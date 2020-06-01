#==============================================================================
#
# Makefile for genomics metric calculation code.
#
#==============================================================================

SHELL   = /bin/sh
#CC      = gcc
CC      = cc

DRIVERDIR = .
SRCDIR = ../src
BINDIR = ../bin
LIBDIR = ../lib

CFLAGS += -Wall -Wno-unused-function -pedantic -std=c99 -Werror
CFLAGS += -I$(SRCDIR) -I$(DRIVERDIR)
CFLAGS += -I../magma/magma_minproduct-1.6.2/include -DADD_
CFLAGS += $(CRAY_CUDATOOLKIT_INCLUDE_OPTS)

CFLAGS_DEBUG = -g
CFLAGS_RELEASE = -DNDEBUG
CFLAGS_OPT = -O2 -ffast-math -fargument-noalias-anything
CFLAGS_OPT += -finline-functions -finline-limit=1000

ifeq ($(FP_PRECISION),single)
  CFLAGS += -DFP_PRECISION_SINGLE
else
  CFLAGS += -DFP_PRECISION_DOUBLE
endif

LFLAGS = $(LDFLAGS)
LFLAGS += -L$(LIBDIR) -lgenomics_metrics
LFLAGS += -L../magma/magma_minproduct-1.6.2/lib -lmagma_minproduct
LFLAGS += $(CRAY_CUDATOOLKIT_POST_LINK_OPTS) -lcublas

EXEC = $(BINDIR)/genomics_metric
LIB = $(LIBDIR)/libgenomics_metrics.a

SOURCES = $(shell echo $(DRIVERDIR)/*.c)
HEADERS = $(shell echo $(DRIVERDIR)/*.h)
OBJECTS = $(SOURCES:.c=.o)

SOURCES_LIB = $(shell echo $(SRCDIR)/*.c)
HEADERS_LIB = $(shell echo $(SRCDIR)/*.h)
OBJECTS_LIB = $(SOURCES_LIB:.c=.o)

.PHONY: debug
debug: CFLAGS += $(CFLAGS_DEBUG)
debug: $(EXEC)

.PHONY: release
release: CFLAGS += $(CFLAGS_RELEASE) $(CFLAGS_OPT)
release: $(EXEC)

#all: $(EXEC)

$(EXEC): $(BINDIR) $(OBJECTS) $(LIB)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LFLAGS)

$(LIB): $(LIBDIR) $(LIB)($(OBJECTS_LIB))
	ranlib $@

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(BINDIR):
	@mkdir -p $(BINDIR)

$(LIBDIR):
	@mkdir -p $(LIBDIR)

clean:
	-rm -f $(OBJECTS)
 
distclean: clean
	-rm -rf $(EXEC) $(BINDIR) $(LIBDIR)

.PHONY: all clean distclean
#==============================================================================
