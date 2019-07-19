CC = g++
CFLAGS = -Wall -O3 --std=c++11

make : main.exe CityGenerator.exe
	
main.exe : main.cpp main.h random.o
	$(CC) main.cpp random.o -o main.exe $(CFLAGS)
	 
CityGenerator.exe : CityGenerator.cpp main.h random.o
	$(CC) CityGenerator.cpp random.o -o CityGenerator.exe $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o *.exe seed.out *.dat
