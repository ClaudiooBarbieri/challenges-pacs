#path to s libraries
ADDITIONAL_PATH := 
EXAMPLES_INCLUDE = /home/claudius/PACS/pacs-examples/Examples/include
MY_INCLUDE = ../include
EXAMPLES_LIB = /home/claudius/PACS/pacs-examples/Examples/lib
# compiler variables

CXXFLAGS = -std=c++20 -O2 -fPIC
CPPFLAGS = -DNDEBUG -I${EXAMPLES_INCLUDE} -I${MY_INCLUDE}
LDFLAGS = -L. -Wl,-rpath=$(ADDITIONAL_PATH) -L${EXAMPLES_LIB}
LDLIBS = 

# files

# source file
SRCS = $(wildcard *.cpp)
# object file from sources
OBJS = $(SRCS:.cpp=.o)
# all headers
HAEDS =$(wildcard *.hpp)
# exe main
EXEC = main

.phony = all clean 

all : $(SRCS) $(OBJS) $(EXEC)

$(EXEC) : $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(EXEC) $(LDLIBS)

$(OBJS) : $(SRCS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(SRCS)

clean :
	$(RM) $(OBJS) $(EXEC)

