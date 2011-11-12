libSigPyProc32.so: libSigPyProc32.o
	ld -shared -L${FFTW_PATH}/lib -L/usr/lib/python2.6/ -lfftw3f -lpython2.6 -soname $@.1 -o $@.1.0 -lc $<
	/sbin/ldconfig -v -n .
	ln -sf $@.1 $@
 
libSigPyProc32.o: SigPyProc32.c Makefile
	gcc -lfftw3f -I/usr/include/python2.6/ -I${FFTW_PATH}/include -L/usr/lib/python2.6/ -L${FFTW_PATH}/lib --fast-math -lm -fPIC -c $< -o $@
 
libSigPyProc8.so: libSigPyProc8.o
	ld -shared -L${FFTW_PATH}/lib -L/usr/lib/python2.6/ -lfftw3f -lpython2.6 -soname $@.1 -o $@.1.0 -lc $< 
	/sbin/ldconfig -v -n .
	ln -sf $@.1 $@
 
libSigPyProc8.o: SigPyProc8.c Makefile
	gcc -lfftw3f -I/usr/include/python2.6/ -I${FFTW_PATH}/include -L/usr/lib/python2.6/ -L${FFTW_PATH}/lib --fast-math -lm -fPIC -c $< -o $@

 
libSigPyProc.so: libSigPyProc.o
	ld -shared -L${FFTW_PATH}/lib -L/usr/lib/python2.6/ -lfftw3f -lpython2.6 -soname $@.1 -o $@.1.0 -lc $< 
	/sbin/ldconfig -v -n .
	ln -sf $@.1 $@
 
libSigPyProc.o: SigPyProc.c Makefile
	gcc -lfftw3f -I/usr/include/python2.6/ -I${FFTW_PATH}/include -L/usr/lib/python2.6/ -L${FFTW_PATH}/lib --fast-math -lm -fPIC -c $< -o $@





all:
	make libSigPyProc32.so
	make libSigPyProc8.so
	make libSigPyProc.so
# Clean target
.PHONY: clean
clean:
	rm -f *.o *.so*
