# 'make'	build executable file
# 'make clean'	removes all *.o and executalbe file

# define the C compiler
CC	= gcc

# define any compile-time flags
CFLAGS = -Wall -g -fPIC -O3 -std=c99

# define openmp flags
OPENMP  = -fopenmp
CUOPENMP  = -Xcompiler -fopenmp

# define the direction containing header file
INCLUDES= -I/usr/local/include -I./

# define the library path
LFLAGS	= -L/usr/local/lib

# define any libraries to link to executable
LIBS	= -lm -lgsl -lgslcblas

# define the C object files
OBJS	=  update.o dtype.o union_find.o models.o stats.o


#define the directory for object
OBJSDIR = object

# define the executable file
MAINSRC = projects/jq3_ladder_square/main.c
MAIN	= exe

all: $(MAIN)

$(MAIN): $(OBJS) main.o
	$(CC) $(CFLAGS) -o $(MAIN) main.o $(OBJS) $(LIBS) $(LFLAGS) $(INCLUDES) 

main.o :
	$(CC) $(CFLAGS) $(INCLUDES) -c $(MAINSRC)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^

lib: $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^



# clean the executable file and object files
clean:
	$(RM) $(OBJS) $(MAIN) main.o
