all: simple gmp

simple:
	$(MAKE) -C src simple
	$(MAKE) -C src clean

gmp:
	$(MAKE) -C src gmp
	$(MAKE) -C src clean

install:
	sudo cp -f gyt-* /usr/local/bin

.PHONY: clean scratch

clean:
	rm -f gyt gyt-*

scratch: clean
	rm -f *~
	rm -f data/*~
