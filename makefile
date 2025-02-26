FILES = bin/cahn-hilliard
CXX = g++
CFLAGS = -w -std=c++17 -I./include -I./include/nlohmann -g # -fsanitize=address 
LDFLAGS = # -fsanitize=address

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:src/%.cpp=obj/%.o)

all: $(FILES)

$(FILES): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

obj/%.o: src/%.cpp | obj
	$(CXX) $(CFLAGS) -c $< -o $@

obj:
	mkdir -p obj

clean:
	rm -f obj/*.o $(FILES)
