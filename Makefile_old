CXXFLAGS = -std=c++11 -fopenmp -O3 -DNDEBUG
CXX = mpic++

INCLUDES =
LDFLAGS =
LIBS =

TARGET = cg
HEAD1 = grid
HEAD2 = solver
OBJS = $(TARGET).o $(HEAD1).o $(HEAD2).o

all: $(TARGET)

$(TARGET): $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)
$(TARGET).o: $(TARGET).cpp Timer.h Makefile
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(TARGET).cpp
$(HEAD1).o: $(HEAD1).h $(HEAD1).cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(HEAD1).cpp
$(HEAD2).o: $(HEAD1).h $(HEAD2).cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(HEAD2).cpp

solve: $(HEAD1).h $(HEAD2).cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(HEAD2).cpp -o $(HEAD2)
	
test:
	mpirun -np 2 ./cg 100 100 2 -1 10 0.1 10 0.5 2

clean:
	@$(RM) -rf *.o *.txt 
