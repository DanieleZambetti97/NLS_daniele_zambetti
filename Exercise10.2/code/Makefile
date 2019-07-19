CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

compila : main.cpp random.cpp CityGenerator.cpp random.h main.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	$(CC) CityGenerator.cpp random.o -o CityGenerator.exe $(CFLAGS)
	./CityGenerator.exe
	$(CC) main.cpp random.o -o main.exe $(CFLAGS)
	mpiexec -np 4 ./main.exe

clean :
	rm *.o *.exe *.dat seed.out 
