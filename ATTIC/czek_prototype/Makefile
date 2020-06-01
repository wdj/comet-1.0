####################################################
# Makefile for Czekanowski Similarity Metric
####################################################
SHELL   = /bin/sh
#CC      = gcc
CC      = cc

#CFLAGS  += -pedantic -Wall -O2
CFLAGS  += -Wall -O2 -finline-functions -finline-limit=1000
CFLAGS  += -I../magma/magma_minproduct-1.6.2/include $(CRAY_CUDATOOLKIT_INCLUDE_OPTS) -DADD_
LFLAGS = $(LDFLAGS) -L../magma/magma_minproduct-1.6.2/lib -lmagma_minproduct $(CRAY_CUDATOOLKIT_POST_LINK_OPTS) -lcublas

ifneq ($(NCOPIES_V),)
  CFLAGS += -DNCOPIES_V=$(NCOPIES_V)
endif

ifneq ($(NCOPIES_F),)
  CFLAGS += -DNCOPIES_F=$(NCOPIES_F)
endif

ifneq ($(PROCESS_VECTORS_FN),)
CFLAGS += -DPROCESS_VECTORS_FN=$(PROCESS_VECTORS_FN)
endif

TARGET  = czek
SOURCES = $(shell echo *.c)
HEADERS = $(shell echo *.h)
OBJECTS = $(SOURCES:.c=.o)

INSTALLDIR  = /usr/local/bin

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.c $(HEADERS)
# %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

install: $(TARGET)
	install -d -m 0755 $(INSTALLDIR)
	install -m 0755 $(TARGET) $(INSTALLDIR)
 
uninstall:
	-rm $(INSTALLDIR)/$(TARGET)

clean:
	-rm -f $(OBJECTS)
 
distclean: clean
	-rm -f $(TARGET)

.PHONY: all install uninstall clean distclean
