FLAGS= --std=c++11 -O3
COMPILER = g++

all: find_MCO_serial.x find_equal_pair.x find_MCO_serial_bulge_tolerant.x #find_MCO_neutral.x 

#find_MCO_neutral.x: find_MCO_neutral.cpp
#	$(COMPILER) $(FLAGS) find_MCO_neutral.cpp -o find_MCO_neutral.x

find_MCO_serial.x: find_MCO_serial.cpp
	$(COMPILER) $(FLAGS) find_MCO_serial.cpp -o find_MCO_serial.x

find_equal_pair.x: find_equal_pair.cpp
	$(COMPILER) $(FLAGS) find_equal_pair.cpp -o find_equal_pair.x
	
find_MCO_serial_bulge_tolerant.x: find_MCO_serial_bulge_tolerant.cpp
	$(COMPILER) $(FLAGS) find_MCO_serial_bulge_tolerant.cpp -o find_MCO_serial_bulge_tolerant.x
	
clean:  
	rm -rf find_MCO_serial.x find_equal_pair.x
 
