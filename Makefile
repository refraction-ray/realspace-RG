CXX = g++
# CXX = icc
MPICC = mpicxx
INC = include
SRC = src
MPISRC = mpisrc
BUILD = build
BIN = bin
TARGET = $(BIN)/rsrg
MPITARGET = $(BIN)/rsrgmpi
CFLAGS = -O2 -I $(INC)
# CFLAGS = -Wall -g -I $(INC) 
# AMARPATH = 
# CFLAGS = -O3 -I $(INC) -I $(AMARPATH) -DARMA_DONT_USE_WRAPPER -mkl=parallel  -openmp  
LIB = -larmadillo 
# LIB = 

SRCTEMP = $(wildcard $(SRC)/*.cpp)
SRCFILE = $(filter-out $(SRC)/main_%.cpp, $(SRCTEMP))
OBJFILE = $(patsubst %.cpp, $(BUILD)/%.o, $(notdir $(SRCFILE)))

MAIN = $(SRC)/main_s.cpp
MPIMAIN = $(SRC)/main_mpi.cpp

# test:
# 	@echo $(SRCFILE)
# 	@echo $(OBJFILE)

all: $(TARGET) $(MPITARGET)

$(MPITARGET): $(OBJFILE) $(MPIMAIN)
	@mkdir -p $(@D)
	$(MPICC) $(CFLAGS) $^ -o $@ $(LIB) -std=c++11

$(TARGET): $(OBJFILE) $(MAIN)
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) $^ -o $@  $(LIB) -std=c++11

$(BUILD)/%.o: $(SRC)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@ -std=c++11

clean:
	rm $(BUILD)/*.o  
