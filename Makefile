BUILD =./build/
BIN =./bin/
CXXFLAGS = -Wall -g -O3
.PHONY = clean, dirs, docs, test1, test2, test3, test4, analyze
CXX = g++

$(BIN)ShallowWater: $(BUILD)ShallowWater.o
	$(CXX) -fopenmp -o $@ $^ -lboost_program_options -lblas

$(BUILD)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -o $@ -c $<

test1:
	$(BIN)ShallowWater --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1 --method 0

test2:
	$(BIN)ShallowWater --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2 --method 1

test3:
	$(BIN)ShallowWater --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --method 1

test4:
	$(BIN)ShallowWater --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4 --method 0

dirs:
	mkdir build bin

clean:
	rm -r build bin html latex

docs:
	doxygen Doxyfile

analyze:
	rm -rf initial.er
	collect -o initial.er bin/ShallowWater --ic 1
	analyzer initial.er