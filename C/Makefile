# 'make'	build executable file
# 'make clean'	removes all *.o and executalbe file

# define the C compiler
CC	= gcc

# define any compile-time flags
CFLAGS = -Wall -g -fPIC -O0 -std=c99

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
OBJS	=  update.o update_omp.o dtype.o union_find.o models.o stats.o


#define the directory for object
OBJSDIR = object

# define the executable file
MAIN	= exe

# define the main srouce file
#MAINSRC = projects/parallel/jq3_ladder_square/main.c
MAINSRC = projects/qlmodel_3d/main.c

all: $(MAIN)

$(MAIN): $(OBJS) main.o
	$(CC) $(CFLAGS) $(OPENMP) -o $(MAIN) main.o $(OBJS) $(LIBS) $(LFLAGS) $(INCLUDES) 

main.o : $(MAINSRC)
	$(CC) $(CFLAGS) $(OPENMP) $(INCLUDES) -c $(MAINSRC)

%.o: %.c
	$(CC) $(CFLAGS) $(OPENMP) $(INCLUDES) -c $^

lib: $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(OPENMP) $(INCLUDES) -c $^



# clean the executable file and object files
clean:
	$(RM) $(OBJS) $(MAIN) main.o
