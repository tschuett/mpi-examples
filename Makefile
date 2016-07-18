CXX=cc
LD=cc
CXXFLAGS=
LDFLAGS=

ifeq ($(VERBOSE),1)
  VC=
  NVC=@true ||
else
  VC=@
  NVC=@
endif

SRC_CPP=$(wildcard *.cpp)

OBJS=$(SRC_CPP:.cpp=.o)

all: cart-pscw.x

%.x: %.o
	$(NVC)echo -e "\e[0;33mCreating $@\033[39m"
	$(LD) $(LDFLAGS) $< -o $@

.c.o: Makefile
	$(NVC)echo -e "\e[0;32mCompiling" $< "\033[39m"
	$(VC)$(CXX) -MD $(CXXFLAGS)  -c $< -o $@

clean:
	@-rm -rf $(shell find . -type f -name '*.o')
	@-rm -rf $(shell find . -type f -name '*.d')
	@-rm -f *.x

-include $(OBJS:.o=.d)

.PHONY: clean
