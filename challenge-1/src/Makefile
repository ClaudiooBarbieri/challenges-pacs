# compiler variables

CXXFLAGS = -std=c++20 -O2 -fPIC
CPPFLAGS = -DNDEBUG -I/home/claudius/PACS/pacs-examples/Examples/include -I../include
LDFLAGS = -L. -Wl,-rpath=${PWD} -L/home/claudius/PACS/pacs-examples/Examples/lib
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

