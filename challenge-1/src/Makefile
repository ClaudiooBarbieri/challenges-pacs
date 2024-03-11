# compiling variables

CXXFLAGS = -std=c++20 -O2 -fPIC
CPPFLAGS = -DNDEBUG -I${mkEigenInc} -I../include
LDFLAGS = -L. -Wl,-rpath=${PWD}
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

