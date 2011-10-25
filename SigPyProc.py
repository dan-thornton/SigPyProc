import os,sys,Header,numpy,pylab
import ctypes as C
import SigPyProcUtils as utils
from scipy.stats.stats import chisquare
lib = C.CDLL("libSigPyProc.so")

#In progress
class BandpassFromBuffer(Header.Header):
    def __init__(self,parentInfo,mallocBuffer=None):
        Header.Header.__init__(self,parentInfo)
        if not mallocBuffer == None:
            self.bpassBuffer = mallocBuffer

class FoldedData(Header.Header):
    def __init__(self,parentInfo,foldBuffer,countsBuffer,
                 period,dm,nbins,nints,nbands):
        Header.Header.__init__(self,parentInfo)
        self.foldBuffer = foldBuffer
        self.countsBuffer = countsBuffer
        self.foldArray = foldBuffer.Ndarray.reshape(nints,nbands,nbins)
        self.countsArray = countsBuffer.Ndarray.reshape(nints,nbands,nbins)
        self.foldArray = self.foldArray/self.countsArray 
        self.period = period
        self.dm = dm
        self.nbins = nbins
        self.nints = nints
        self.nbands = nbands
        self.reset()
        
    def _relevel(self):
        bad_ids = numpy.where(numpy.isnan(self.foldArray))
        good_ids = numpy.where(self.countsArray>0)
        self.foldArray[bad_ids] = numpy.median(self.foldArray[good_ids])
        self.reset()
        
    def reset(self):
        self.profile = self.foldArray.sum(axis=0).sum(axis=0)
        self.freqPhase = self.foldArray.sum(axis=0)
        self.timePhase = self.foldArray.sum(axis=1)

    def getDMdelays(self,dm):
        chanwidth = self.info['foff']*self.info["nchans"]/self.nbands
        chanFreqs = (numpy.arange(self.nbands)
                     *chanwidth)+self.info['fch1']
        chanDelays = (dm * 4.148808e3 *
                      ((chanFreqs**-2)-(self.info['fch1']**-2)
                       )/(self.info["tsamp"]/self.period)
                      ).round().astype("int32")%self.nbins
        return chanDelays
        
    def getPdelays(self,p):
        Wsubint = self.nints/(self.info["tsamp"]*self.info["nsamples"])
        Wbin = self.period/self.nbins
        drifts = p*(numpy.arange(self.nints)
                    *Wsubint/Wbin).round()%self.nbins
        return drifts.astype("int32")

    def optimiseDM(self,hidm=50,lodm=-50,dmstep=1):
        dms = numpy.arange(lodm,hidm,float(dmstep))
        dmcurve = numpy.empty_like(dms)
        for jj,dm in enumerate(dms):
            nFreqPhase = numpy.empty_like(self.freqPhase)
            delays = self.getDMdelays(dm)
            for ii in range(self.nbands):
                nFreqPhase[ii] = utils.rollArray(self.freqPhase[ii],delays[ii],0)
            dmcurve[jj] = chisquare(nFreqPhase.sum(axis=0))[0]
        return dmcurve

    def optimiseP(self,loP=-0.005,hiP=0.005,Pstep=0.000005):
        periods = numpy.arange(loP,hiP,float(Pstep))
        pcurve = numpy.empty_like(periods)
        for jj,p in enumerate(periods):
            nTimePhase = numpy.empty_like(self.timePhase)
            delays = self.getPdelays(p)
            for ii in range(self.nints):
                nTimePhase[ii] = utils.rollArray(self.timePhase[ii],delays[ii],0)
            pcurve[jj] = chisquare(nTimePhase.sum(axis=0))[0]
        return pcurve

    def plot(self,cmap="jet",interp="gaussian"):
        pylab.subplot(221)
        pylab.title("Time Vs Phase")
        pylab.xlabel("Phase")
        pylab.ylabel("Subintegration")
        pylab.imshow(self.timePhase,aspect="auto",cmap=cmap,interpolation=interp)
        
        pylab.subplot(222)
        pylab.title("Frequency Vs Phase")
        pylab.xlabel("Phase")
        pylab.ylabel("Subband")
        pylab.imshow(self.freqPhase/self.freqPhase.sum(axis=1).reshape(self.nbands,1),
                     aspect="auto",cmap=cmap,interpolation=interp)

        pylab.subplot(223)
        pylab.xlabel("Phase")
        pylab.ylabel("Power")
        pylab.plot(self.profile/self.profile.mean())
       
        if self.nbands > 1:
            chiDM = self.optimiseDM()
            pylab.subplot(426)
            pylab.xlabel("DM")
            pylab.ylabel("SNR")
            pylab.plot(chiDM)

        if self.nints > 1:
            chiP = self.optimiseP()
            pylab.subplot(428)
            pylab.xlabel("Period")
            pylab.ylabel("SNR")
            pylab.plot(chiP)

        pylab.show()
        


class TimeSeries(Header.Header):
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
      
        Header.Header.__init__(self,parentInfo)
        if isinstance(parentInfo,dict) and not mallocBuffer == None:
            self.timBuffer = mallocBuffer
    
        elif isinstance(parentInfo,str) and mallocBuffer == None:
            self._readTim()
    

    def _readTim(self):
        """Read time series from .tim file.
        
        Currently only works with 32-bit time series.
        Creates class attribute timBuffer.
        """

        f = utils.File(self.file,"r")
        f.seek(self.hdrlen)
        self.timBuffer = utils.Buffer(self.info["nsamples"],C.c_float)
        f.cread(self.timBuffer)

    def fold(self,period,nbins=50,nints=1):
        """Fold time series into phase bins and subintegrations.

        Args:
        period -- period in seconds to fold with
        nbins  -- number of phase bins in output
        nsubs  -- number of subintegrations in output
        
        Returns: FoldedData instance  
        """

        foldBuffer  = utils.Buffer(nbins*nints,C.c_double)
        countBuffer = utils.Buffer(nbins*nints,C.c_int)
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
        fftBuffer = utils.Buffer(fftsize+2,C.c_float)
        lib.rfft(self.timBuffer.Cbuffer,fftBuffer.Cbuffer,fftsize)
        return SpectraFromBuffer(self.info.copy(),fftBuffer)

    def runningMean(self,window=10001):
        """Filter time series by a running mean.

        Args:
        window -- width in bins of running mean filter
        """
        newTimBuffer  = utils.Buffer(self.info["nsamples"],C.c_float)
        lib.runningMean(self.timBuffer.Cbuffer,newTimBuffer.Cbuffer,self.info["nsamples"],window)
        return TimeSeries(self.info.copy(),newTimBuffer)
        

    def toPrestoFormat(self,basename):
        """Write time series in presto .dat format.
        """
        self.make_inf(outfile="%s.inf"%(basename))
        datfile = utils.File("%s.dat"%(basename),"w+")
        self.timBuffer.Ndarray.toFile(datfile)
        
    

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
        print minlag 
        for key in keys:
            self.lags[key] += abs(minlag)
             
        


class SpectraFromBuffer(Header.Header):
    """Class to handle output of fftw.

    Parameters
    ----------

    \tparentInfo   -- Dictionary containing header info 
    \tmallocBuffer -- SigPyProcUtils.Buffer instance containing complex spectrum

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
        Header.Header.__init__(self,parentInfo)
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
        timBuffer = utils.Buffer(timsize,C.c_float)
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
        allCands = utils.stackRecarrays(tuple(allCands))

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
                    

class Filterbank(Header.Header):
    """Class to handle filterbank files in all their glory.

    Parameters
    ----------

    \tfilename   -- name of filterbank file

    Methods
    -------

    \tcollapse     -- Collapse filterbank in frequency and/or time 
                      (ie. DM0 time series and/or bandpass)
    \tdedisperse   -- Simple dedispersion algorithm (fixed some sigproc bugs)
    \tdownsample   -- Downsample filterbank in frequency and/or time
    """

    def __init__(self,filename):
        """Create Filterbank instance.

        Args:
        filename -- string containing name of file to process
        """
        Header.Header.__init__(self,filename)
        self.f = utils.File(filename,"r")

    def getDMdelays(self,dm):
        chanFreqs = (numpy.arange(self.info["nchans"])
                     *self.info['foff'])+self.info['fch1']

        chanDelays = (dm * 4.148808e3 *
                      ((chanFreqs**-2)-(self.info['fch1']**-2)
                       )/self.info["tsamp"]).round().astype("int32")

        return chanDelays

    def collapse(self,blocksize=512,freq=True,time=True):
        """Get bandpass and/or zero-dm timeseries from data.

        Args:
        blocksize -- Number of samples to read in for each block.
        freq      -- Flag True to collapse in frequency
        time      -- Flag True to collapse in time

        Returns: 
        if time and freq: tuple of BandpassFromBuffer and TimeSeries instances
        if time and not freq: TimeSeries instance
        if freq and not time: BandpassFromBuffer instance
        """

        self.f.seek(self.hdrlen)
        readBuffer = utils.Buffer(blocksize*self.info["nchans"],C.c_ubyte)
        if time: timBuffer = utils.Buffer(self.info["nsamples"],C.c_float)
        if freq: bpassBuffer = utils.Buffer(self.info["nchans"],C.c_float)
        lastread = self.info["nsamples"]%blocksize
        nreads = (self.info["nsamples"]-lastread)/blocksize
        for ii in xrange(nreads):
            self.f.cread(readBuffer,blocksize*self.info["nchans"])
            if time: lib.getTim(readBuffer.Cbuffer,timBuffer.Cbuffer,
                                self.info["nchans"],blocksize,ii*blocksize)
            if freq: lib.getBpass(readBuffer.Cbuffer,bpassBuffer.Cbuffer,
                                  self.info["nchans"],blocksize)

        self.f.cread(readBuffer,lastread*self.info["nchans"])
        if time: lib.getTim(readBuffer.Cbuffer,timBuffer.Cbuffer,
                            self.info["nchans"],lastread,(ii+1)*blocksize)
        if freq: lib.getBpass(readBuffer.Cbuffer,bpassBuffer.Cbuffer,
                              self.info["nchans"],lastread)

        if time and freq: return(BandpassFromBuffer(self.info.copy(),bpassBuffer),
                                 TimeSeries(self.info.copy(),timBuffer))
        elif time and not freq: return TimeSeries(self.info.copy(),timBuffer)
        elif freq and not time: return BandpassFromBuffer(self.info.copy(),bpassBuffer)
        else: return None

    def dedisperse(self,dm,gulp=10000):
        """Dedisperse filterbank to timeseries.

        Args:
        dm        -- Dispersion measure to dedisperse to 
        gulp      -- size of block to read at a time, if chosen gulp
                     is less than maximum dispersion delay gulp is taken
                     as 2 * max delay.

        Returns: TimeSeries instance
        """

        chanDelays = self.getDMdelays(dm)
        delayPointer = utils.arrayAsPointer(chanDelays)
        maxDelay = int(chanDelays.max())
        gulp = max(2*maxDelay,gulp)
        timLen = self.info["nsamples"]-maxDelay

        self.f.seek(self.hdrlen)
        timBuffer = utils.Buffer(timLen,C.c_float)
        readBuffer = utils.Buffer(self.info["nchans"]*gulp,C.c_ubyte)
        nreads = int(self.info["nsamples"]/(gulp-maxDelay))
        lastread = self.info["nsamples"]%(gulp-maxDelay)
        
        for ii in xrange(nreads):
            self.f.cread(readBuffer,gulp*self.info["nchans"])
            lib.dedisperse(readBuffer.Cbuffer,timBuffer.Cbuffer,delayPointer,
                           maxDelay, self.info["nchans"], gulp, ii*(gulp-maxDelay))
            self.f.seek(self.hdrlen+(self.info["nchans"]*(ii+1)*(gulp-maxDelay)))
        self.f.cread(readBuffer,lastread*self.info["nchans"])
        lib.dedisperse(readBuffer.Cbuffer,timBuffer.Cbuffer,delayPointer,
                       maxDelay, self.info["nchans"], lastread, (ii+1)*(gulp-maxDelay))
        timInfo = self.info.copy()
        timInfo["nsamples"] = timLen
        timInfo["refdm"] = dm
        return TimeSeries(timInfo,timBuffer)

    
    def downsample(self,tfactor=1,ffactor=1,filename="Outfile.fil"):
        """Downsample filterbank in frequency and/or time.

        Args:
        tfactor  -- Factor to downsample in time
        ffactor  -- Factor to downsample in frequecy (must be factor of nchans)
        filename -- File to write downsampled data to.

        Returns: Filterbank instance (from output file)
        """

        self.f.seek(self.hdrlen)
        if not self.info["nchans"]%ffactor == 0:
            print "Bad frequency factor given"
            return None

        self.downsample_header(tfactor=tfactor,ffactor=ffactor)
        outFile = utils.File(filename,"w+")
        outFile.write(self.write_header())
        self.reset_header()
        
        readBuffer = utils.Buffer(tfactor*self.info["nchans"],C.c_ubyte)
        writeBuffer = utils.Buffer(self.info["nchans"]/ffactor,C.c_ubyte)
        tempBuffer = utils.Buffer(self.info["nchans"]/ffactor,C.c_int)

        lib.downsampleFil(self.f.f, outFile.f, readBuffer.Cbuffer,
                          writeBuffer.Cbuffer,tempBuffer.Cbuffer,
                          tfactor, ffactor, self.info["nchans"],
                          self.info["nsamples"])

        return Filterbank(filename)

    def fold(self,period,dm,nbins=50,nints=32,nbands=32,gulp=10000):
        """Fold filterbank phase bins, subintegrations and subbands. 

        Args:
        period -- period in seconds to fold with
        nbins  -- number of phase bins in output
        nints  -- number of subintegrations in output
        
        

        Returns: FoldedData instance
        """

        nbands = min(nbands,self.info["nchans"])
        chanDelays = self.getDMdelays(dm) 
        delayPointer = utils.arrayToPointer(chanDelays)
        maxDelay = int(chanDelays.max())
        gulp = max(2*maxDelay,gulp)

        self.f.seek(self.hdrlen)
        
        foldBuffer  = utils.Buffer(nbins*nints*nbands,C.c_float)
        countBuffer = utils.Buffer(nbins*nints*nbands,C.c_int)
        readBuffer = utils.Buffer(self.info["nchans"]*gulp,C.c_ubyte)
        
        nreads = int(self.info["nsamples"]/(gulp-maxDelay))
        lastread = self.info["nsamples"]%(gulp-maxDelay)
        
        for ii in xrange(nreads):
            sys.stdout.write("%d%%\r"%(100*ii/nreads))
            sys.stdout.flush()
            self.f.cread(readBuffer,gulp*self.info["nchans"])
            lib.foldFil(readBuffer.Cbuffer, foldBuffer.Cbuffer, countBuffer.Cbuffer,
                        delayPointer, maxDelay, C.c_double(self.info["tsamp"]),
                        C.c_double(period), gulp, self.info["nchans"], nbins,
                        nints, nbands, ii*(gulp-maxDelay))

            self.f.seek(self.hdrlen+(self.info["nchans"]*(ii+1)*(gulp-maxDelay)))

        self.f.cread(readBuffer,lastread*self.info["nchans"])

        lib.foldFil(readBuffer.Cbuffer, foldBuffer.Cbuffer, countBuffer.Cbuffer,
                    delayPointer, maxDelay, C.c_double(self.info["tsamp"]), 
                    C.c_double(period), gulp, self.info["nchans"], nbins,
                    nints, nbands, (ii+1)*(gulp-maxDelay))

        parentInfo = self.info.copy()
        return FoldedData(parentInfo,foldBuffer,countBuffer,
                          period,dm,nbins,nints,nbands)
        

class MultiFilterbank:
    def __init__(self,files):
        self.files = files
        self.filterbanks = [Filterbank(filename) for filename in files]
        
            
        
