import numpy,pylab,sys,cPickle
import ctypes as C
from Header import Header
from SigPyProcUtils import File,Buffer,arrayToPointer
lib = C.CDLL("libSigPyProc.so")


class Filterbank(Header):
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
        Header.__init__(self,filename)
        self.f = File(filename,"r")

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

        if self.info["nbits"] == 32:
            readBuffer = Buffer(blocksize*self.info["nchans"],C.c_float)
            getTim = lib.getTim32
            getBpass = lib.getBpass32
        elif self.info["nbits"] == 8:
            readBuffer = Buffer(blocksize*self.info["nchans"],C.c_ubyte)
            getTim = lib.getTim8
            getBpass = lib.getBpass8

        self.f.seek(self.hdrlen)
        if time: timBuffer = Buffer(self.info["nsamples"],C.c_float)
        if freq: bpassBuffer = Buffer(self.info["nchans"],C.c_float)
        lastread = self.info["nsamples"]%blocksize
        nreads = (self.info["nsamples"]-lastread)/blocksize
        for ii in xrange(nreads):
            sys.stdout.write("%d%%\r"%(100*ii/nreads))
            sys.stdout.flush()
            self.f.cread(readBuffer,blocksize*self.info["nchans"])
            if time: getTim(readBuffer.Cbuffer,timBuffer.Cbuffer,
                                self.info["nchans"],blocksize,ii*blocksize)
            if freq: getBpass(readBuffer.Cbuffer,bpassBuffer.Cbuffer,
                                  self.info["nchans"],blocksize)

        self.f.cread(readBuffer,lastread*self.info["nchans"])
        if time: getTim(readBuffer.Cbuffer,timBuffer.Cbuffer,
                            self.info["nchans"],lastread,(ii+1)*blocksize)
        if freq: getBpass(readBuffer.Cbuffer,bpassBuffer.Cbuffer,
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
        delayPointer = arrayToPointer(chanDelays)
        maxDelay = int(chanDelays.max())
        gulp = max(2*maxDelay,gulp)
        timLen = self.info["nsamples"]-maxDelay

        self.f.seek(self.hdrlen)
        timBuffer = Buffer(timLen,C.c_float)

        if self.info["nbits"] == 32:
            readBuffer = Buffer(self.info["nchans"]*gulp,C.c_float)
            dedisperser = lib.dedisperse32
        elif self.info["nbits"] == 8:
            readBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
            dedisperser = lib.dedisperse8

        nreads = int(self.info["nsamples"]/(gulp-maxDelay))
        lastread = self.info["nsamples"]%(gulp-maxDelay)
        
        for ii in xrange(nreads):
            sys.stdout.write("%d%%\r"%(100*ii/nreads))
            sys.stdout.flush()
            self.f.cread(readBuffer,gulp*self.info["nchans"])
            dedisperser(readBuffer.Cbuffer,timBuffer.Cbuffer,delayPointer,
                           maxDelay, self.info["nchans"], gulp, ii*(gulp-maxDelay))
            self.f.seek(self.hdrlen+(self.info["nchans"]*(ii+1)*(gulp-maxDelay)))
        self.f.cread(readBuffer,lastread*self.info["nchans"])
        dedisperser(readBuffer.Cbuffer,timBuffer.Cbuffer,delayPointer,
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
        outFile = File(filename,"w+")
        outFile.write(self.write_header())
        self.reset_header()
        
        readBuffer = Buffer(tfactor*self.info["nchans"],C.c_ubyte)
        writeBuffer = Buffer(self.info["nchans"]/ffactor,C.c_ubyte)
        tempBuffer = Buffer(self.info["nchans"]/ffactor,C.c_int)

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
        delayPointer = arrayToPointer(chanDelays)
        maxDelay = int(chanDelays.max())
        gulp = max(2*maxDelay,gulp)

        self.f.seek(self.hdrlen)
        
        foldBuffer  = Buffer(nbins*nints*nbands,C.c_float)
        countBuffer = Buffer(nbins*nints*nbands,C.c_int)
        readBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
        
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
        
    
    def chanStats(self,window=10001):        

        self.f.seek(self.hdrlen)
        gulp = 10*window
        lastread = self.info["nsamples"]%(gulp-window)
        nreads = (self.info["nsamples"]-lastread)/(gulp-window)
        

        if self.info["nbits"] == 32:
            readBuffer = Buffer(self.info["nchans"]*gulp,C.c_float)
            writeBuffer = Buffer(self.info["nchans"]*gulp,C.c_float)
            maximaBuffer = Buffer(self.info["nchans"],C.c_float)
            minimaBuffer = Buffer(self.info["nchans"],C.c_float)
            getStats = lib.getStats32
        elif self.info["nbits"] == 8:
            readBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
            writeBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
            maximaBuffer = Buffer(self.info["nchans"],C.c_ubyte)
            minimaBuffer = Buffer(self.info["nchans"],C.c_ubyte)
            getStats = lib.getStats8

        meansBuffer = Buffer(self.info["nchans"],C.c_double)
        bpassBuffer = Buffer(self.info["nchans"],C.c_double)
        stdevBuffer = Buffer(self.info["nchans"],C.c_double)

        outFile = File("%s_RM.fil"%(self.basename),"w+")
        outFile.write(self.write_header())
        infoFile = File("%s_RM.info"%(self.basename),"w+")
        
        for ii in xrange(nreads):
            sys.stdout.write("%d%%\r"%(100*ii/nreads))
            sys.stdout.flush()
            self.f.cread(readBuffer,gulp*self.info["nchans"])
            
            getStats(readBuffer.Cbuffer, meansBuffer.Cbuffer, bpassBuffer.Cbuffer, stdevBuffer.Cbuffer,
                     writeBuffer.Cbuffer,maximaBuffer.Cbuffer,minimaBuffer.Cbuffer,self.info["nchans"],gulp,window,ii)
            
            if ii == 0:
                outFile.cwrite(writeBuffer)
            else:
                outFile.cwrite(writeBuffer, nunits=(gulp-window)*self.info["nchans"],
                               offset=self.info["nchans"]*window)

            self.f.seek(self.hdrlen+4*self.info["nchans"]*((ii+1)*(gulp-window)))
            
        self.f.cread(readBuffer,lastread*self.info["nchans"])
        getStats(readBuffer.Cbuffer, meansBuffer.Cbuffer, bpassBuffer.Cbuffer, stdevBuffer.Cbuffer,
                 writeBuffer.Cbuffer,maximaBuffer.Cbuffer,minimaBuffer.Cbuffer,self.info["nchans"],lastread,window,ii)

        outFile.cwrite(writeBuffer, nunits=(lastread-window)*self.info["nchans"],
                       offset=self.info["nchans"]*window)
        stdevBuffer.Ndarray = numpy.sqrt(stdevBuffer.Ndarray/self.info["nsamples"])
        
        info = {"sigma":stdevBuffer.Ndarray,
                "bandpass":bpassBuffer.Ndarray,
                "maxima":maximaBuffer.Ndarray,
                "minima":minimaBuffer.Ndarray}
        cPickle.dump(info,infoFile)

        
class MultiFilterbank:
    def __init__(self,files):
        self.files = files
        self.filterbanks = [Filterbank(filename) for filename in files]




        
from TimeSeries import TimeSeries
from FoldedData import FoldedData
from Bandpass import BandpassFromBuffer            
        
