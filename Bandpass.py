import numpy,pylab
import ctypes as C
from SigPyProcUtils import Buffer,File
from Header import Header
from scipy.signal import medfilt
lib = C.CDLL("libSigPyProc.so")

#In progress
class BandpassFromBuffer(Header):
    def __init__(self,parentInfo,mallocBuffer=None):
        Header.__init__(self,parentInfo)
        if not mallocBuffer == None:
            self.bpassBuffer = mallocBuffer

