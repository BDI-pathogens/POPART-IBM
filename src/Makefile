#	compiling instructions:
#	make clean
#	make all
# 
#	For compilation on Rescomp:
#	make clean
#	make all location=rescomp
#	
#	For compilation on Rescomp using the icc compiler:
#	make clean
#	make all location=rescomp compiler=icc


# Executing instructions:
#	./popart-simul.exe Zambia

GSLFLAGS= -lgsl -lgslcblas -lm -O3
LFLAGS = ${GSLFLAGS}

# all the C files (*.c)
SRCS = constants.c demographics.c hiv.c init.c input.c interventions.c output.c partnership.c utilities.c checks.c memory.c fitting.c debug.c pc.c simul.c main.c 

# all the object files (*.o)
OBJS = constants.o demographics.o hiv.o init.o input.o interventions.o output.o partnership.o utilities.o checks.o memory.o fitting.o debug.o pc.o simul.o main.o 

# the executable (*.exe)
_EXE = popart-simul.exe 
EXE = $(_EXE)

# For use of the icc compiler
ifeq ($(compiler),icc)
        C = icc
else
        C = gcc
endif

# For compilation on Rescomp
ifeq ($(location),rescomp)
	INC = /apps/well/gsl/2.2.1-gcc4.9.3/include
	LIB = /apps/well//gsl/2.2.1-gcc4.9.3/lib
else
        INC = /usr/local/include
	LIB = /usr/local/lib
endif

# Compilation options and libraries to be used
CFLAGS = -g -fcommon -Wall -fmessage-length=0 -I$(INC) -O0

# To compile
all: ${OBJS}
	${C} -L${LIB} -o ${EXE} ${OBJS} ${LFLAGS}

clean:
	rm -f $(OBJS) $(EXE)

.c.o:
	${C} ${CFLAGS} -c $< -o $@	
