# Driver code
DRIVER = generateMIPCuts
EPENG = EquitablePartition
AGGENG = Aggregate
ORBREAD = read_orbits

EXE = $(DRIVER)
COINDIR=/home/edeakins/miniconda3

SRCDIR=/home/edeakins/MIPSymmetryCuts/CutGeneration
# the compiler: gcc for C program, define as g++ for C++
CXX = g++ --std=c++11
CXXLINKFLAGS =  -Wl,--rpath -Wl,/home/edeakins/miniconda3/lib
CXXFLAGS = -g -pedantic-errors

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings

INC = -I$(COINDIR)/include/coin -I$(SRCDIR)/
LIB = -L$(COINDIR)/lib -lClpSolver -lClp -lCoinUtils -lCgl -lOsi -lOsiClp -lz -lm 

# The build target 

all: $(EXE)

$(EXE): $(DRIVER).cpp
	$(CXX) $(CXXFLAGS) $(CXXLINKFLAGS) $(INC) -o $(EXE) $(DRIVER).cpp $(EPENG).cpp $(AGGENG).cpp $(ORBREAD).cpp $(LIB) $(ADDLIBS)

clean:
	rm -rf generateMIPCuts

# $(TARGET): $(TARGET).c
#     $(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c