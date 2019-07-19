CC = g++
CFLAGS = -Wall -O3 --std=c++11
main.exe : datablocking.o MolDyn_NVE.o
	$(CC) MolDyn_NVE.o datablocking.o -o MolDyn_NVE.x
MolDyn_NVE.o: MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
datablocking.o : datablockingspecial.cpp datablockingspecial.h
	$(CC) -c datablockingspecial.cpp -o datablocking.o $(CFLAGS)
clean: 
	rm -rf risultati/* frames/*.xyz *.o *.x
esegui:
	./MolDyn_NVE.x
