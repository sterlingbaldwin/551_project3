CC = icpc
FLAGS = -g -Wall -fopenmp -O3 -o
TARGET = gauss
SRC = gauss.cpp

gauss: $(SRC)
	$(CC) $(FLAGS) $(TARGET) $(SRC)
