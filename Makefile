all: SVF

SVF:
	make -C src -f Makefile1 install

clean:
	make -C src -f Makefile1 clean
	rm -f bin/*.exe
