
# Variable for path to PACS Examples
PACS_ROOT = /home/claudius/PACS/pacs-examples/Examples

# Path for reaching libraries and header providede by pacs-Examples repository
EXAMPLES_INCLUDE = ${PACS_ROOT}/include
EXAMPLES_LIB = ${PACS_ROOT}/lib

# Path for personal header file
MY_INCLUDE = ../include

# Adding run time search path
ADDITIONAL_PATH = 

# compiler variables

CXXFLAGS = -std=c++20 -O3 -fPIC
CPPFLAGS = -DNDEBUG -I${EXAMPLES_INCLUDE} -I${MY_INCLUDE}
LDFLAGS = -L. -Wl,-rpath=$(EXAMPLES_LIB) -L${EXAMPLES_LIB}
LDLIBS = 

# files

# source file
SRCS = $(wildcard *.cpp)
# object file from sources
OBJS = $(SRCS:.cpp=.o)
# all headers
HEADS =$(wildcard ${MY_INCLUDE}/*.hpp)
# exe main
EXEC = main

.phony = all clean 

# compiling routine
all : $(SRCS) $(OBJS) $(EXEC) 

$(EXEC) : $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(EXEC) $(LDLIBS)

$(OBJS) : $(SRCS) $(HEADS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(SRCS)

# removing obj file and main
clean :
	$(RM) $(OBJS) $(EXEC)
