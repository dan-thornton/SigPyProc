import os,sys,numpy
import ctypes as C
from HeaderParams import nptypes_to_ctypes as np_to_c
clib = C.CDLL("libc.so.6")

class Malloc:
    """Class to wrap c function malloc.
    
    Parameters
    ----------
    \tnunits -- number of units of size sizeof(dtype) in buffer
    \tdtype  -- ctypes type designation
    \tsetbuf  -- flag for performing memset(0) in constructor

    Methods
    -------
    \tsetBuffer -- memset buffer to a given value
    \ttoNdarray -- return the buffer as a numpy array instance 
                   !!!!(be very very careful with reference counts here)!!!!
    \ttoFile    -- write buffer to a Cfile instance
    """

    def __init__(self,nunits,dtype,setbuf=True):
        """Create new Malloc instance by calling malloc from c library.

        Args:
        nunits -- number of units of size sizeof(dtype) in buffer
        dtype  -- ctypes type designation
        setbuf  -- flag for performing memset to zero (this is default behaviour)
        """

        self.dtype = dtype
        self.nunits = int(nunits)
        self.nbytes = C.sizeof(self.dtype)*self.nunits
        clib.malloc.restype = C.POINTER(self.dtype)
        self.buffer = clib.malloc(self.nbytes)
        if setbuf:
            self.setBuffer(0)

    def setBuffer(self,value):
        """Set buffer to a value.

        Args:
        value -- value (of type dtype) to set buffer to
        """

        clib.memset(self.buffer,value,self.nbytes)

    def toNdarray(self):
        """Return malloc'd buffer as Ndarray.

        Warning: this directly accesses the malloc'd buffer as a numpy array, meaning
        that if the block is free'd before you have finished with the array, the array
        contents will be unallocated memory, leading to possible segfaults.
        """

        Carray = self.dtype * self.nunits
        Narray = Carray.from_address(C.addressof(self.buffer.contents))
        return numpy.ctypeslib.as_array(Narray)

    def toFile(self,CfileInst,nbytes=None):
        """Write nbytes of buffer to a Cfile instance.

        Args:
        CfileInst -- A Cfile instance
        nbytes    -- number of bytes to write (default = whole buffer)
        """

        if not nbytes:
            nbytes = self.nbytes
        CfileInst.write(self.buffer,nbytes)

    def __del__(self):
        """Clean up the buffer.

        Uses the free method when the reference count of the malloc buffer reaches 0.
        """
        
        print "Freeing %d * %s byte buffer"%(self.nunits,self.dtype.__name__)
        clib.free(self.buffer)


class Cfile:
    """Wrapper for c FILE* types.

    Wrapper is designed to provide the python file interface, but to c files.
    
    Parameters
    ----------
    \tfilename -- name of file to open
    \tmode     -- mode of access 

    Methods
    -------
    \tseek  -- performs fseek
    \tread  -- performs fread
    \twrite -- performs fwrite
    \tell   -- performs ftell
    """    
    
    def __init__(self,filename,mode):
        """Create new Cfile instance

        Args:
        filename -- name of file to open                                                                                                                               
        mode     -- mode of access (can be any normal c file mode, eg. "r","r+","a" etc.)
        """
        self.filename = filename
        self.f = clib.fopen(self.filename,mode)
        clib.fread.restype = C.c_int

    def __repr__(self):
        return "C file object <%s>"%(self.filename)

    def seek(self,bytes,pos=0):
        """Seek to position in file.

        Args:
        bytes -- bytes to move to from pos
        pos   -- point to move from (pos=0=SEEK_START, pos=1=SEEK_CURRENT, pos=2=SEEK_END)
        """
        clib.fseek(self.f,int(bytes),pos)

    def read(self,membuffer,nbytes):
        """Read data into a given malloc buffer.

        Args:
        membuffer -- A Malloc instance buffer attribute
        nbytes    -- Number of bytes to read
        """
        return clib.fread(membuffer,1,int(nbytes),self.f)

    def write(self,membuffer,nbytes):
        """write data from a given buffer.
        
        Args:
        membuffer -- A Malloc instance buffer attribute                                                                                                                                                           
        nbytes    -- Number of bytes to write
        """
        clib.fwrite(membuffer,1,int(nbytes),self.f)

    def tell(self):
        """Find file position.
        """
        return clib.ftell(self.f)

    def __del__(self):
        """Close file when reference count drops to zero.
        """
        print "Closing file %s"%(self.filename)
        clib.fclose(self.f)

def rollArray(y,shift,axis):
    """Roll the elements in the array by 'shift' positions along the                                                                                                                                                         
    given axis.
    
    Args:
    y       -- array to roll
    shift   -- number of bins to shift by
    axis    -- axis to roll along

    Returns: shifted Ndarray 
    """
    y = numpy.asanyarray(y)
    n = y.shape[axis]
    shift %= n 
    return y.take(numpy.concatenate((numpy.arange(shift,n),numpy.arange(shift))), axis)

def stackRecarrays(arrays):
    """Wrapper for stacking numpy.recarrays.
    
    Args:
    arguments are an arbitrary number of numpy.recarrays with the same dtpyes argument in a tuple
    
    Note: This is required as numpy.hstack will not stack recarrays
    """
    return arrays[0].__array_wrap__(numpy.hstack(arrays))

def nearestFactor(n,val):
    """Find nearest factor.
    
    Args:
    n     -- number that we wish to factor
    val   -- number that we wish to find nearest factor to
    """

    fact=[1,n]
    check=2
    rootn=numpy.sqrt(n)
    while check<rootn:
        if n%check==0:
            fact.append(check)
            fact.append(n/check)
        check+=1
    if rootn==check:
        fact.append(check)
    fact.sort()
    return fact[numpy.abs(numpy.array(fact)-val).argmin()]

def arrayToPointer(array):
    """Return numpy array as a pointer.
    
    Args:
    array -- numpy.array instance
    """
    if not np_to_c.get(array.dtype.str) == None:
        return array.ctypes.data_as(C.POINTER(np_to_c[array.dtype.str]))
    else:
        raise KeyError,"Key '%s' not found in ndtypes_to_ctypes conversion dictionary in HeaderParams module"

