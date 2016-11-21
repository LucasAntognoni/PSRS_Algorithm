SRC = $(wildcard src/*.c)
OBJ = $(patsubst src/%.c, obj/%.o, $(SRC))

CC = mpicc
CFLAGS = -std=c11 -Wall -Wextra -fopenmp -g 
LDFLAGS = -lm

obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(LIBS) -c $< -o $@

all: $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $(OBJ) -o psrs 