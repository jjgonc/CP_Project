CC = gcc
BIN = bin/
SRC = src/
INCLUDES = include/
EXEC = k_means
CFLAGS = -O2 -fopenmp -ftree-vectorize -mavx
.DEFAULT_GOAL = k_means
THREADS = 40
k_means: $(SRC)k_means.c $(BIN)utils.o
	$(CC)  $(SRC)k_means.c $(CFLAGS) $(BIN)utils.o -o $(BIN)$(EXEC)

$(BIN)utils.o: $(SRC)utils.c $(INCLUDES)utils.h
	$(CC) $(CFLAGS) -c $(SRC)utils.c -o $(BIN)utils.o
clean:
	rm -r bin/*
runseq:
	./$(BIN)$(EXEC) 10000000 $(CP_CLUSTERS)
runpar:
	./$(BIN)$(EXEC) 10000000 $(CP_CLUSTERS) $(THREADS)