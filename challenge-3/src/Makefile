# Variable for path to PACS Examples
PACS_ROOT = /home/claudius/PACS/pacs-examples/Examples

# Path for reaching libraries and header provided by pacs-Examples repository
EXAMPLES_INCLUDE = ${PACS_ROOT}/include
EXAMPLES_LIB = ${PACS_ROOT}/lib

# Path for personal header file
MY_INCLUDE = ../include

# Adding run time search path
ADDITIONAL_PATH = 

# Compiler variables
CXXFLAGS = -fopenmp -std=c++20 -O3 -fPIC 
CPPFLAGS = -DNDEBUG -I${EXAMPLES_INCLUDE} -I${MY_INCLUDE}
LDFLAGS = -L. -Wl,-rpath=$(EXAMPLES_LIB) -L${EXAMPLES_LIB}
LDLIBS = -lmuparser -lgomp

# MPI compiler variables
MPI_CXX = mpic++ -fopenmp
MPI_LDLIBS = -lmuparser -lmpi -lgomp

# Files
# Source file
SRCS = $(wildcard *.cpp)
# Object files from sources
OBJS = $(SRCS:.cpp=.o)
# All headers
HEADS = $(wildcard ${MY_INCLUDE}/*.hpp)
# Executable names
EXEC = main
MPI_EXEC = main_mpi

.PHONY: all clean mpi

# Compiling routine
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(EXEC) $(LDLIBS) $(MPI_LDLIBS)

# Pattern rule for object files
%.o: %.cpp $(HEADS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# parallel compiling routine
mpi: $(MPI_EXEC)

$(MPI_EXEC): $(OBJS)
	$(MPI_CXX) $(LDFLAGS) $(OBJS) -o $(MPI_EXEC) $(LDLIBS) $(MPI_LDLIBS)

clean:
	$(RM) $(OBJS) $(EXEC) $(MPI_EXEC)
