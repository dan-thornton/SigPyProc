import numpy,pylab
import ctypes as C
from Header import Header
#from TimeSeries import TimeSeries
from SigPyProcUtils import File,Buffer,stackRecarrays
lib = C.CDLL("libSigPyProc.so")

class SpectraFromBuffer(Header):
    """Class to handle output of fftw.

    Parameters
    ----------

    \tparentInfo   -- Dictionary containing header info 
    \tmallocBuffer -- SigPyProcBuffer instance containing complex spectrum

    Methods
    -------

    \tiFFt         -- Perform inverse complex to real 1D DFT (uses fftw3)
    \trednoiseSinc -- Multipy red end of spectrum by sum of sinc functions (experimental)
    \tsumHarmonics -- Sum harmonics through downsampling and summing
    \tfindCands    -- Find all peaks above a given sigma value in the spectrum
    """

    def __init__(self,parentInfo,mallocBuffer):
        """Construct spectrum object from buffer.

        Args:
        parentInfo -- dict instance containing header information
        mallocBuffer -- a SigPyProcUtils Buffer instance.

        """
        Header.__init__(self,parentInfo)
        self.fftBuffer = mallocBuffer
        self.info = parentInfo
        self.candDtype = [("bin","int"),
                          ("freq","float"),
                          ("period","float"),
                          ("sigma","float"),
                          ("fold","int"),
                          ("refdm","float")]
        self._formSpec()
        
    def _formSpec(self):
        """Form power spectrum from real and complex fft output.
        """
        ar = self.fftBuffer.Ndarray.astype("float64").\
            reshape((self.info["nsamples"]+2)/2,2).\
            transpose()
        self.realPart  = ar[0]
        self.imagPart  = ar[1]
        self.specArray = ar[0]**2 + ar[1]**2
        self.harmonicSums = None

    def _getSinc(self,ii,nbins=100):
        """Get a sinc function.
        
        Args:
        ii    -- bin of first minimum
        nbins -- length of output
 
        Returns: 1-sinc fuction of length nbins
        """

        maxVal = (nbins/float(ii))*numpy.pi/2.
        stepSize = numpy.pi/(ii*2.)
        func = 1-numpy.sinc(numpy.arange(0,maxVal,stepSize))
        if func.shape[0]%nbins == 0:
            return func
        else:
            return func[:-1]

    def iFFT(self):
        """Perform 1D complex to real inverse FFT using fftw3.

        Returns: TimeSeries instance
        """

        timsize = (self.specArray.shape[0]*2)-2
        timBuffer = Buffer(timsize,C.c_float)
        lib.ifft(self.fftBuffer.Cbuffer,timBuffer.Cbuffer,timsize)
        return TimeSeries(self.info.copy(),timBuffer)
        
    def rednoiseSinc(self,nbins=100):
        """Multiply part of spectrum by sum of 1-sinc functions.

        Args:
        nbins -- spectrum[0:nbins] will be multiplied by the sum of 
                 sinc functions with first minima between 2 and nbins
        """

        sincCurve = numpy.array([self._getSinc(ii,nbins) for ii in xrange(2,nbins)])\
            .sum(axis=0)/(nbins-2)
        self.specArray[0:nbins] = self.specArray[0:nbins]*sincCurve
        
    def sumHarmonics(self,numharms=4):
        """Perform harmonic folding via a squeezing routine.

        Args:
        numharms -- number of times to fold (def=4, ie. 16th harmonic)
        """

        self.harmonicSums = []
        array = self.specArray
        for ii in xrange(numharms):
            array = self._squeezeSum(array)
            self.harmonicSums.append(array)
                    
    def _squeezeSum(self,array):
        """Squeeze array and sum.
        
        Args:
        array -- array to sum

        Returns: numpy.array instance
        """

        if not array.shape[0]%2 == 0: array = array[:-1]
        return array[:array.shape[0]/2]+ array.reshape(array.shape[0]/2,2).sum(axis=1)

    def findCands(self,sigma=6):
        """Find candidate signals in power spectrum and harmonic folds.

        Args:
        sigma -- minimum sigma value to return
        
        Returns: numpy.recarray instance (sorted (by sigma) list of cands)
        """

        allCands = []
        allCands.append(self._getPeaks(sigma,self.specArray,0))

        if not self.harmonicSums == None:
            for ii,spectrum in enumerate(self.harmonicSums):
                allCands.append(self._getPeaks(sigma,spectrum,ii+1))
        allCands = stackRecarrays(tuple(allCands))

        sortedCands = numpy.sort(allCands,order="sigma")[::-1]

        print "Period(s)  Freq(Hz)  Sigma  Bin  Fold"
        for line in sortedCands:
            print "%g\t  %g\t  %g\t  %g\t  %g"%\
                (line["period"],line["freq"],line["sigma"],
                 line["bin"],line["fold"])
        return sortedCands 
            
    def _getPeaks(self,sigma,spectrum,fold=0):
        """Find peaks in power spectrum. 

        Args:
        sigma    -- minimum sigma value to return as peak
        spectrum -- spectrum to search for signals
        fold     -- fold number (required for correct period/frequency determination)
        
        Returns: numpy.recarray instance (candidates)
        """

        if spectrum==None: spectrum = self.specArray 
        sigmaArray = (spectrum-spectrum.mean())/spectrum.std()
        idx = numpy.where(sigmaArray>sigma)
        nCands = idx[0].shape[0]
        print "Found %d cands above %f sigma"%(nCands,sigma)
        candArray = numpy.recarray(nCands,dtype=self.candDtype)
        candArray["bin"] = idx[0]
        candArray["sigma"] = sigmaArray[idx]
        candArray["freq"] = (idx[0]*1/(self.info["tsamp"]*self.info["nsamples"]))/(2**fold)
        candArray["period"] = 1/candArray["freq"]
        candArray["fold"] = fold
        candArray["refdm"] = self.info.get("refdm",0.0)
        return candArray
                    
from TimeSeries import TimeSeries
