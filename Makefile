CXX = g++
LDFLAGS = -larmadillo

OBJ = $(SRC:src/%.cpp=build/%.o)
SRC = $(wildcard src/*.cpp)
EXE = programa

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

build/%.o: src/%.cpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf build/*

