SRC = $(wildcard src/*.c)
WORKER_SRC = $(wildcard src/worker/*.c)

CC = mpicc
CFLAGS = -std=c11 -g -Wall -Wextra 
LDFLAGS = -lm

all:
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $(SRC) -fopenmp -o psrs
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $(WORKER_SRC) -o worker 