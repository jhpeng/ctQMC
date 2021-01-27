# 'make'	build executable file
# 'make clean'	removes all *.o and executalbe file

# define the C compiler
CC	= gcc

# define any compile-time flags
CFLAGS = -Wall -Werror -Wextra -g -O0 -std=c99

# define openmp flags
OPENMP  = -fopenmp
CUOPENMP  = -Xcompiler -fopenmp

# define the direction containing header file
INCLUDES= -I/usr/local/include -I./

# define the library path
LFLAGS	= -L/usr/local/lib 

# define any libraries to link to executable
LIBS	= -lm -lgsl -lgslcblas 


#define the directory for object
OBJS = dtype.o union_find.o models.o update.o stats.o

# define the executable file
MAIN = lib/libctqmc.a

all: $(MAIN)

$(MAIN): $(OBJS)
	ar -rcs $(MAIN) *.o 
	cp *.h include/

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^

# clean the executable file and object files
clean:
	$(RM) $(OBJS) $(MAIN)
