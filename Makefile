BUILD =./build/
BIN =./bin/
CXXFLAGS = -Wall -g
.PHONY = clean, dirs, docs, test1, test2, test3, test4
CXX = g++

$(BIN)ShallowWater: $(BUILD)ShallowWater.o
	$(CXX) -o $@ $^ -lboost_program_options

$(BUILD)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

test1:
	$(BIN)ShallowWater --dT 0.1 --T 80 --Nx 100 --Ny 100 --ic 1

test2:
	$(BIN)ShallowWater --dT 0.1 --T 80 --Nx 100 --Ny 100 --ic 2

test3:
	$(BIN)ShallowWater --dT 0.1 --T 80 --Nx 100 --Ny 100 --ic 3

test4:
	$(BIN)ShallowWater --dT 0.1 --T 80 --Nx 100 --Ny 100 --ic 4

dirs:
	mkdir build bin

clean:
	rm -r build bin html latex

docs:
	doxygen Doxyfile