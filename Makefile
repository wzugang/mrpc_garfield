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

gas_properties: gas_properties.C 
	$(CXX) $(CFLAGS) gas_properties.C
	$(CXX) -o gas_properties gas_properties.o $(LDFLAGS)
	rm gas_properties.o

