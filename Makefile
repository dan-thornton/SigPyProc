# This is a minial makefile for building a shared library
 
# test.so is a shared library built from test.c
libSigPyProc.so: libSigPyProc.o
	ld -shared -L${FFTW_PATH}/lib -lfftw3f -soname $@.1 -o $@.1.0 -lc $<
	/sbin/ldconfig -v -n .
	ln -sf $@.1 $@
 
# This is an intermediate object. Must be built with
# the -fPIC option
libSigPyProc.o: SigPyProc.c Makefile
	gcc -lfftw3f -I${FFTW_PATH}/include -L${FFTW_PATH}/lib -lm -fPIC -c $< -o $@
 
# Clean target
.PHONY: clean
clean:
	rm -f *.o

