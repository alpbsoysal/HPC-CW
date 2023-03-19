BUILD =./build/
BIN =./bin/
CXXFLAGS = -Wall
.PHONY = clean, dirs, run
CXX = g++

$(BIN)ShallowWater: $(BUILD)ShallowWater.o
	$(CXX) -o $@ $^ -lboost_program_options

$(BUILD)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

run:
	$(BIN)ShallowWater

dirs:
	mkdir build bin

clean:
	rm build/* bin/*