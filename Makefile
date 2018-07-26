CXXFLAGS += -g -Wall -O3 --std=c++11

OBJS = mergeSNPlogs diploidizeSNPlog compareSNPlogs

.PHONY: all,clean

all: mergeSNPlogs diploidizeSNPlog compareSNPlogs

clean:
	rm $(OBJS)
