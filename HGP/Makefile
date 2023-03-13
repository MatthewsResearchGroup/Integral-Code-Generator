CXX=g++
CXXFLAGS=-std=c++17 -I../common


exe:=hrr
src:=$(shell find . -iname "*cxx") $(shell find ../common -iname "*.cxx")
objs:=$(addsuffix .o, $(basename $(src)))
deps:=$(join $(addsuffix .deps/,$(dir $(objs))),$(notdir $(objs:.o=.d)))

.PHONY: all
all: $(exe)

-include $(deps)

.PHONY: clean
clean:
	rm -f $(objs) $(exe)

$(exe): $(objs)
	@mkdir -p $(dir $(exe))
	$(CXX) $(LDFLAGS) -o $(exe) $(objs) 

%.o: %.cxx
	@mkdir -p $(dir $@).deps
	$(CXX) $(CXXFLAGS) -MT $@ -MF $(dir $@).deps/$(notdir $(@:.o=.d)) -o $@ -c $<


