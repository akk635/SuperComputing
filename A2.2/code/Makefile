####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

CC = mpicc
CFLAGS = -Wall -g -O0 -std=c99
LIBS = -lm -lmetis
METIS_INCLUDE = $(METIS_INC) 
METIS_LIB_PATH = $(METIS_LIB)
CFLAGS += $(PAPI_INC)
LIBS += $(PAPI_LIB)


LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

SRCS = initialization.c compute_solution.c finalization.c test_functions.c util_read_files.c util_write_files.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: gccg 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(METIS_INCLUDE) $(INCLUDE)

gccg: gccg.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) $(METIS_INCLUDE) $(METIS_LIB_PATH) $(LIBS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o gccg $(LIBPOS) *.vtk *.dSYM
