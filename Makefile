OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library
GFORTRANDIR = /usr/local/Cellar/gfortran/4.8.2/gfortran/lib

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR)
#	-O3 \

# Debug flags
CFLAGS += -g

LDFLAGS = `root-config --glibs` -lGeom -L$(GFORTRANDIR) -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield
LDFLAGS += -g

mrpc: mrpc.C 
	$(CXX) $(CFLAGS) $<
	$(CXX) -o $@ $@.o $(LDFLAGS)
	rm $@.o

gasfile: gasfile.C 
	$(CXX) $(CFLAGS) gasfile.C
	$(CXX) -o gasfile gasfile.o $(LDFLAGS)
	rm gasfile.o

