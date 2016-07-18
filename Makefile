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

CFLAGS=

all: cart-pscw.x

%.x: %.o
	$(NVC)echo -e "\e[0;33mCreating $@\033[39m"
	$(VC)$(LD) $(LDFLAGS) $< -o $@

.c.o: Makefile
	$(NVC)echo -e "\e[0;32mCompiling" $< "\033[39m"
	$(VC)$(CC) $(CFLAGS) -c $< -o $@

clean:
	@-rm -rf $(shell find . -type f -name '*.o')
	@-rm -f *.x
