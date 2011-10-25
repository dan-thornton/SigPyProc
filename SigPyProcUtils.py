import os,sys,numpy
import ctypes as C
from HeaderParams import nptypes_to_ctypes as np_to_c
clib = C.CDLL("libc.so.6")
lib  = C.CDLL("libSigPyProc.so") 

class Buffer:
    """Class to wrap c array allocation.
    
    Parameters
    ----------
    \tnunits -- number of units of size sizeof(dtype) in buffer
    \tdtype  -- ctypes type designation

    Methods
    -------
    \tsetBuffer -- memset buffer to a given value
    """

    def __init__(self,nunits,dtype):
        """Create new Buffer instance.
        
        Args:
        nunits -- number of units of size sizeof(dtype) in buffer
        dtype  -- ctypes type designation
        """
        self.dtype = dtype
        self.nunits = int(nunits)
        self.nbytes = C.sizeof(self.dtype)*self.nunits
        plan = self.nunits * self.dtype
        self.Cbuffer = plan()
        self.Ndarray = numpy.ctypeslib.as_array(self.Cbuffer) 

    def setBuffer(self,value):
        """Set buffer to a value.

        Args:
        value -- value (of type dtype) to set buffer to
        """
        clib.memset(self.Cbuffer,value,self.nbytes)

    def __del__(self):
        """Clean up the buffer.
        
        Uses the free method when the reference count of the malloc buffer reaches 0.
        """
        print "Freeing %d * %s byte buffer"%(self.nunits,self.dtype.__name__)


class File(file):
    """
    Wrapper is designed to provide the python file interface, but to c files.
    
    Parameters
    ----------
    \tfilename -- name of file to open
    \tmode     -- mode of access 

    Methods
    -------
    \tall the methods of file 
    """    
    
    def __init__(self,filename,mode):
        """Create new File instance

        Args:
        filename -- name of file to open
        mode     -- mode of access 
        """
        file.__init__(self,filename,mode)
        lib.cfile.argtypes = [C.py_object]
        self.cfile = lib.cfile(self)
        self.filename = filename

    def cread(self,BufInst,nunits=None):
        if nunits == None:
            clib.fread(BufInst.Cbuffer,C.sizeof(BufInst.dtype),BufInst.nunits,self.cfile)
        else:
            clib.fread(BufInst.Cbuffer,C.sizeof(BufInst.dtype),nunits,self.cfile)
        
    def __del__(self):
        """Close file when reference count drops to zero.
        """
        print "Closing file %s"%(self.filename)
        self.close()


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

