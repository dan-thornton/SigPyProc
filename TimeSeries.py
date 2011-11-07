import numpy,pylab
import ctypes as C
from Header import Header
from SigPyProcUtils import File,Buffer
lib = C.CDLL("libSigPyProc.so")

class TimeSeries(Header):
    """Class to handle time series data.

    Parameters
    ----------
    
    \tparentInfo   -- Dictionary containing header info or 32bit .tim file
    \tmallocBuffer -- SigPyProcUtils.Buffer instance containing time series 
    
    Methods
    -------
    
    \tfold           -- Fold a time series
    \trFFt           -- Perform forward real to complex 1-d DFT (uses fftw3)  
    \trunningMean    -- Remove any inherent baseline in the data
    \ttoPrestoFormat -- Write time series into presto compatible format (with inf file)
    """

    def __init__(self,parentInfo,mallocBuffer=None):
        """Construct time series object from buffer or file.

        Args:
        parentInfo -- may be a file name or a dict instance containing 
                      header information

        mallocBuffer -- a SigPyProcUtils Buffer instance.
                        If parentInfo is a file, this argument
                        is ignored
        """ 
      
        Header.__init__(self,parentInfo)
        if isinstance(parentInfo,dict) and not mallocBuffer == None:
            self.timBuffer = mallocBuffer
    
        elif isinstance(parentInfo,str) and mallocBuffer == None:
            self._readTim()
            

    def _readTim(self):
        """Read time series from .tim file.
        
        Currently only works with 32-bit time series.
        Creates class attribute timBuffer.
        """

        f = File(self.file,"r")
        f.seek(self.hdrlen)
        self.timBuffer = Buffer(self.info["nsamples"],C.c_float)
        f.cread(self.timBuffer)

    def fold(self,period,nbins=50,nints=1):
        """Fold time series into phase bins and subintegrations.

        Args:
        period -- period in seconds to fold with
        nbins  -- number of phase bins in output
        nsubs  -- number of subintegrations in output
        
        Returns: FoldedData instance  
        """

        foldBuffer  = Buffer(nbins*nints,C.c_double)
        countBuffer = Buffer(nbins*nints,C.c_int)
        lib.foldTim(self.timBuffer.Cbuffer,foldBuffer.Cbuffer,
                    countBuffer.Cbuffer,C.c_double(self.info["tsamp"]),
                    C.c_double(period),self.info["nsamples"],nbins,nints)
        parentInfo = self.info.copy()
        return FoldedData(parentInfo,foldBuffer,countBuffer,
                          period,self.info.get("refdm",0),
                          nbins,nints,1)

    def rFFT(self):
        """Perform 1D real to complex forward FFT using fftw3.

        Returns: SpectraFromBuffer instance  
        """

        if self.info["nsamples"]%2 ==0:
            fftsize = self.info["nsamples"]
        else:
            fftsize = self.info["nsamples"]-1
        fftBuffer = Buffer(fftsize+2,C.c_float)
        lib.rfft(self.timBuffer.Cbuffer,fftBuffer.Cbuffer,fftsize)
        return SpectraFromBuffer(self.info.copy(),fftBuffer)

    def runningMean(self,window=10001):
        """Filter time series by a running mean.

        Args:
        window -- width in bins of running mean filter
        """
        newTimBuffer  = Buffer(self.info["nsamples"],C.c_float)
        lib.runningMean(self.timBuffer.Cbuffer,newTimBuffer.Cbuffer,self.info["nsamples"],window)
        return TimeSeries(self.info.copy(),newTimBuffer)
        

    def toPrestoFormat(self,basename):
        """Write time series in presto .dat format.
        """
        self.make_inf(outfile="%s.inf"%(basename))
        datfile = File("%s.dat"%(basename),"w+")
        self.timBuffer.Ndarray.tofile(datfile)
        
    

class MultiTimeSeries:
    def __init__(self,files=None):
        self.tims = {}
        self.lags = {}
        if isinstance(files,list):
            for filename in files:
                self.tims[filename] = TimeSeries(filename)
                                
    def append(self,TimeSeriesInst):
        self.tims[TimeSeriesInst.info["filename"]] = TimeSeriesInst
    
    def findLags(self,window=200):
        keys = self.tims.keys()
        reference = self.tims[keys[0]]
        self.lags[keys[0]] = 0
        for key in keys[1:]:
            corr = numpy.correlate(reference.timBuffer.Cbuffer[window:-window],
                                   self.tims[key].timBuffer.Cbuffer,mode="valid")
            self.lags[key] = corr.argmax()-window
        minlag = min(self.lags.values())
        for key in keys:
            self.lags[key] += abs(minlag)
             
from FoldedData import FoldedData
from Spectra import SpectraFromBuffer
