OBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp) $(wildcard em/*.cpp))

include ../Mk/local.Mk
CXXINCLUDE+=-I./ -I../common 
CXXINCLUDE+=-I/usr/local/include

all: $(OBJS) 

.cpp.o:
	$(CXX)  $(CXXFLAGS) $(CXXDEFINE) -c  $(CXXINCLUDE) $< 

.c.o:
	$(CXX)  $(FLAGS) $(CXXDEFINE) -c  $(CXXFLAGS) $(CXXINCLUDE) $< 

clean:
	$(RM) $(OBJS) *.moc $(UI_HEAD) $(UI_CPP)
