CC = g++

CCFLAGS = -O3 -static
OBJS = mdl.o

all	: mdl

clean   :
	rm *.o mdl

mdl  :	${OBJS}
	${CC} -o $@ ${CCFLAGS} ${OBJS}

mdl.o :	mdl.h mdl.c
	${CC} -c ${CCFLAGS} mdl.c

