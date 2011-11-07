import os,numpy,pylab,sys,cPickle,threading,time
import ctypes as C
from Header import Header,MultiHeader
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
        infofile = "%s.info"%(self.basename)
        if os.path.isfile("%s.info"%(self.basename)):
            try:
                self.stats = cPickle.load(open(infofile,"r"))
            except:
                self.stats = None
        else:
            self.stats = None
            


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

        nreads = self.info["nsamples"]//(gulp-maxDelay)
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
        
        nreads = self.info["nsamples"]//(gulp-maxDelay)
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
        
    
    def getStatistics(self,window=10001):        

        self.f.seek(self.hdrlen)
        gulp = 2*window
        lastread = self.info["nsamples"]%(gulp-window)
        nreads = self.info["nsamples"]//(gulp-window)
        

        if self.info["nbits"] == 32:
            readBuffer = Buffer(self.info["nchans"]*gulp,C.c_float)
            writeBuffer = Buffer(self.info["nchans"]*gulp,C.c_float)
            maximaBuffer = Buffer(self.info["nchans"],C.c_float)
            minimaBuffer = Buffer(self.info["nchans"],C.c_float)
            getStats = lib.getStats32
        elif self.info["nbits"] == 8:
            raise ValueError,"8bit mode currently not working"
            #readBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
            #writeBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
            #maximaBuffer = Buffer(self.info["nchans"],C.c_float)
            #minimaBuffer = Buffer(self.info["nchans"],C.c_float)
            #getStats = lib.getStats8

        meansBuffer = Buffer(self.info["nchans"],C.c_double)
        bpassBuffer = Buffer(self.info["nchans"],C.c_double)
        stdevBuffer = Buffer(self.info["nchans"],C.c_double)

        outFile = File("%s_RM.fil"%(self.basename),"w+")
        outFile.write(self.write_header())
        
        for ii in xrange(nreads):
            sys.stdout.write("%d%%\r"%(100*ii/nreads))
            sys.stdout.flush()
            self.f.cread(readBuffer,gulp*self.info["nchans"])
            
            getStats(readBuffer.Cbuffer, meansBuffer.Cbuffer,
                     bpassBuffer.Cbuffer, stdevBuffer.Cbuffer,
                     writeBuffer.Cbuffer,maximaBuffer.Cbuffer,
                     minimaBuffer.Cbuffer,self.info["nchans"],gulp,window,ii)
            
            if ii == 0:
                outFile.cwrite(writeBuffer)
            else:
                outFile.cwrite(writeBuffer, nunits=(gulp-window)*self.info["nchans"],
                               offset=self.info["nchans"]*window)

            self.f.seek(self.hdrlen+4*self.info["nchans"]*((ii+1)*(gulp-window)))
            
        self.f.cread(readBuffer,lastread*self.info["nchans"])
        getStats(readBuffer.Cbuffer, meansBuffer.Cbuffer, 
                 bpassBuffer.Cbuffer, stdevBuffer.Cbuffer,
                 writeBuffer.Cbuffer,maximaBuffer.Cbuffer,
                 minimaBuffer.Cbuffer,self.info["nchans"],lastread,window,ii)

        outFile.cwrite(writeBuffer, nunits=(lastread-window)*self.info["nchans"],
                       offset=self.info["nchans"]*window)
        stdevBuffer.Ndarray = numpy.sqrt(stdevBuffer.Ndarray/self.info["nsamples"])
        
        infoFile = open("%s_RM.info"%(self.basename),"w+")
        info = {"sigma":stdevBuffer.Ndarray,
                "bandpass":bpassBuffer.Ndarray,
                "maxima":maximaBuffer.Ndarray,
                "minima":minimaBuffer.Ndarray}
        cPickle.dump(info,infoFile)
        infoFile.close()
        return Filterbank(outFile.name)
        
    def dropBits32to8(self,gulp=1024,flag=0,clip=0):
        if self.stats is None:
            raise AttributeError,"Filterbank instance has no attribute 'stats'"
        
        lastread = self.info["nsamples"]%gulp
        nreads = self.info["nsamples"]//gulp
        
        readBuffer = Buffer(self.info["nchans"]*gulp,C.c_float)
        writeBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
        
        if flag is not 0:
            flagBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
            flagMax = flag*self.stats["sigma"].astype("float32")
            flagMin = -flag*self.stats["sigma"].astype("float32")
            
        if clip is not 0:
            chanMax = numpy.array([min(a,b) for a,b in zip(clip*self.stats["sigma"],self.stats["maxima"])])
            chanMin = numpy.array([max(a,b) for a,b in zip(-clip*self.stats["sigma"],self.stats["minima"])])
        else:
            chanMax = self.stats["maxima"]
            chanMin = self.stats["minima"]

        chanFactor = ((chanMax-chanMin)/255.).astype("float32")
        chanPlus = (chanMin/chanFactor).astype("float32")
        
        outFile = File("%s_8bit.fil"%(self.basename),"w+")
        self.modify_header("nbits",8)
        outFile.write(self.write_header())
        if flag is not 0:
            flagFile = File("%s_8bit_mask.fil"%(self.basename),"w+")
            self.modify_header("source_name","%s_mask"%(self.info["source_name"].split()[0]))
            flagFile.write(self.write_header())        
        self.reset_header()
        
        self.f.seek(self.hdrlen)

        for ii in xrange(nreads):
            sys.stdout.write("%d%%\r"%(100*ii/nreads))
            sys.stdout.flush()
            self.f.cread(readBuffer,gulp*self.info["nchans"])
            lib.dropBits32to8(readBuffer.Cbuffer,writeBuffer.Cbuffer,
                              arrayToPointer(chanFactor),arrayToPointer(chanPlus),
                              gulp,self.info["nchans"])
            outFile.cwrite(writeBuffer)

            if flag is not 0: 
                lib.flagBlock(readBuffer.Cbuffer,flagBuffer.Cbuffer,
                              arrayToPointer(flagMax),arrayToPointer(flagMin),
                              gulp,self.info["nchans"])
                flagFile.cwrite(flagBuffer)
            
        self.f.cread(readBuffer,lastread*self.info["nchans"])

        lib.dropBits32to8(readBuffer.Cbuffer,writeBuffer.Cbuffer,
                          arrayToPointer(chanFactor),arrayToPointer(chanPlus),
                          lastread,self.info["nchans"])
        outFile.cwrite(writeBuffer,nunits=lastread*self.info["nchans"])

        if flag is not 0:
            lib.flagBlock(readBuffer.Cbuffer,flagBuffer.Cbuffer,
                          arrayToPointer(flagMax),arrayToPointer(flagMin),
                          lastread,self.info["nchans"])
            flagFile.cwrite(flagBuffer,nunits=lastread*self.info["nchans"])

        if flag is not 0:
            return(Filterbank(outFile.name),Filterbank(flagFile.name))
        else:
            return Filterbank(outFile.name)

              
class MultiFilterbank(MultiHeader):
    def __init__(self,files):
        MultiHeader.__init__(self,files)
        self.filterbanks = [Filterbank(filename) for filename in files]
        self.nfiles = len(files)
        self.lags = numpy.zeros(self.nfiles)
        
    def findLags(self):
        zeroDMs = MultiTimeSeries()
        for fil in self.filterbanks:
            zeroDMs.append(fil.collapse(freq=False))
        zeroDMs.findLags()
        self.lags = numpy.array([zeroDMs.lags[fil.info["filename"]] for fil in self.filterbanks])
        return zeroDMs,self.lags

    def genMBmask(self,gulp=512,threshold=4):
        outFile = File("MBmask.fil","w+")
        clags = Buffer(self.lags.shape[0],C.c_int)
        clags.Ndarray = self.lags
        for fil,offset in zip(self.filterbanks,self.lags):
            fil.f.seek(fil.hdrlen+offset)
        minlen = min([ fil.info["nsamples"]-offset for fil,offset in zip(self.filterbanks,self.lags)])
        
        readBuffers = Buffer(gulp*self.info["nchans"],C.c_ubyte,dim=self.nfiles)
        writeBuffer = Buffer(gulp*self.info["nchans"],C.c_ubyte)

        lastread = minlen%gulp
        nreads = minlen//gulp
        
        for ii in xrange(nreads):
            for jj,fil in enumerate(self.filterbanks):
                fil.f.cread(readBuffers,n=jj)
                
            lib.genMBmask(readBuffers.Cbuffer,writeBuffer.Cbuffer,
                          threshold,self.nfiles,gulp,self.info["nchans"])
            outFile.cwrite(writeBuffer)

        for ii,fil in enumerate(self.filterbanks):
            fil.f.cread(readBuffers,lastread*self.info["nchans"],n=ii)

        lib.genMBmask(readBuffers.Cbuffer,writeBuffer.Cbuffer,
                      threshold,self.nfiles,lastread,self.info["nchans"])
        outFile.cwrite(writeBuffer,lastread*self.info["nchans"])
                
from TimeSeries import TimeSeries,MultiTimeSeries
from FoldedData import FoldedData
from Bandpass import BandpassFromBuffer            
        
