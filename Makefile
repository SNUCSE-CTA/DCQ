# CXX := g++
CXX := ~/gcc/bin/g++

CXXFLAGS := -std=c++0x -O3 -w #-DDEBUG
CPPFLAGS := -Iinclude -D_METHOD_C #-DBREAKDOWN #-DVERBOSE

SRC := src
OBJ := obj

SOURCES := $(wildcard $(SRC)/*.cpp)
OBJECTS := $(SOURCES:$(SRC)/%.cpp=$(OBJ)/%.o)

TARGET := DCQ

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.cpp | $(OBJ)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ):
	mkdir -p $(OBJ)

clean:
	@$(RM) -rv $(TARGET) $(OBJ)

