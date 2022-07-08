#FLAGS = -std=c++11 -DARMA_DONT_USE_WRAPPER -larmadillo -llapack -lblas -openmp -O3 -pg

O_FLAGS= -std=c++11 -DARMA_DONT_USE_WRAPPER -I/${ARMA_INC} -O3  -openmp
L_FLAGS= -std=c++11 -DARMA_DONT_USE_WRAPPER -L/${ARMA_LIB} -O3 -llapack -lblas  -openmp

COMPILER = mpic++ #/opt/homebrew/Cellar/open-mpi/4.1.2/bin/mpic++

evolution: evolution.o functions.o random.o
	$(COMPILER) $(L_FLAGS) evolution.o functions.o random.o -o evolution
 
evolution.o: evolution.cpp functions.cpp functions.h random.cpp random.h
	$(COMPILER) -c $(O_FLAGS) evolution.cpp -o evolution.o
	
functions.o: functions.cpp functions.h random.cpp random.h
	$(COMPILER) -c $(O_FLAGS) functions.cpp -o functions.o
	
random.o: random.cpp random.h
	$(COMPILER)  -c $(O_FLAGS) random.cpp -o random.o

clean:  
	rm -r -f evolution histogram* pred* selection *.o
 
