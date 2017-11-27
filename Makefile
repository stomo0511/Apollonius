CXX =	/usr/local/bin/g++ -fopenmp
# for DEBUG
#CXXFLAGS =	-DDEBUG -g
# for Performance evaluation
CXXFLAGS =	-O3

OBJS =		Apollonius.o

all:ap

ap: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

clean:
	rm -f VC $(OBJS)
