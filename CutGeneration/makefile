# Driver code
DRIVER = generateMIPCuts
EPENG = EquitablePartition
AGGENG = Aggregate

EXE = $(DRIVER)

SRCDIR = .
# the compiler: gcc for C program, define as g++ for C++
CXX = g++ --std=c++11
CXXLINKFLAGS =  -Wl,--rpath -Wl,/home/edeakins/Cbc-2.9.8/lib
CXXFLAGS = -g -O0 -pipe -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DCLP_BUILD

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings

INC = -I/home/edeakins/Cbc-2.9.8/include/coin -I/home/edeakins/EQLPSolver/MIPCuts/
LIB = -L/home/edeakins/Cbc-2.9.8/lib -lClpSolver -lClp  -L/home/edeakins/Cbc-2.9.8/lib -lCoinUtils -L/home/edeakins/Cbc-2.9.8/lib -lCgl -L/home/edeakins/Cbc-2.9.8/lib -lOsi -L/home/edeakins/Cbc-2.9.8/lib -lOsiClp -lz -lm 

# The build target 

all: $(EXE)

$(EXE): $(DRIVER).cpp
	$(CXX) $(CXXFLAGS) $(CXXLINKFLAGS) $(INC) -o $(EXE) $(DRIVER).cpp $(EPENG).cpp $(AGGENG).cpp $(LIB) $(ADDLIBS)

clean:
	rm -rf generateMIPCuts

# $(TARGET): $(TARGET).c
#     $(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c