CC = g++
CFLAGS = -Wall -O3 --std=c++11


main.exe : main.o random.o metropolis.o datablocking.o
	$(CC) random.o metropolis.o datablocking.o main.cpp -o main.exe
main.o: main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
datablocking.o : datablocking.cpp datablocking.h
	$(CC) -c datablocking.cpp -o datablocking.o $(CFLAGS)

 
isto: main_isto.o random.o metropolis.o datablocking.o
	$(CC) random.o metropolis.o datablocking.o main_isto.cpp -o main_isto.exe
main_isto: main_isto.cpp
	$(CC) -c main_isto.cpp -o main_isto.o $(CFLAGS)	

clean :
	rm *.o main.exe main_isto.exe seed.out Risultati.dat Posizioni.dat Energie.dat
