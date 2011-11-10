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
        """For a given dispersion measure return delay for each channel
        """
        chanFreqs = (numpy.arange(self.info["nchans"])
                     *self.info['foff'])+self.info['fch1']
        chanDelays = (dm * 4.148808e3 *
                      ((chanFreqs**-2)-(self.info['fch1']**-2)
                       )/self.info["tsamp"]).round().astype("int32")
        return chanDelays

    def collapse(self,blocksize=512):
        """Collapse a filterbank in frequency.

        Args:
        blocksize -- number of samples to read in each gulp
        
        Return: TimeSeries instance
        """
        readBuffer = Buffer(blocksize*self.info["nchans"],C.c_ubyte)
        timBuffer = Buffer(self.info["nsamples"],C.c_float)
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readOnce(blocksize)
        for nsamps,ii in passPlan.makePass():
            lib.getTim8(readBuffer.Cbuffer,timBuffer.Cbuffer,
                        self.info["nchans"],nsamps,ii*blocksize)
        return TimeSeries(self.info.copy(),timBuffer)

    def bandpass(self,blocksize=512):
        """Collapse a filterbank in time.

        Args:
        blocksize -- number of samples to read in each gulp
        
        Return: BandpassFromBuffer instance
        """
        readBuffer = Buffer(blocksize*self.info["nchans"],C.c_ubyte)
        bpassBuffer = Buffer(self.info["nchans"],C.c_float)
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readOnce(blocksize)
        for nsamps,ii in passPlan.makePass():
            lib.getBpass8(readBuffer.Cbuffer,bpassBuffer.Cbuffer,
                          self.info["nchans"],blocksize)
        return BandpassFromBuffer(self.info.copy(),bpassBuffer)

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
        timBuffer = Buffer(timLen,C.c_float)
        readBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readSkipBack(gulp,maxDelay)
        for nsamps,ii in passPlan.makePass():
            lib.dedisperse8(readBuffer.Cbuffer,timBuffer.Cbuffer,delayPointer,
                            maxDelay, self.info["nchans"], nsamps, ii*(gulp-maxDelay))
        timInfo = self.info.copy()
        timInfo["nsamples"] = timLen
        timInfo["refdm"] = dm
        return TimeSeries(timInfo,timBuffer)
    
    def downsample(self,tfactor=1,ffactor=1,filename=None):
        """Downsample filterbank in frequency and/or time.

        Args:
        tfactor  -- Factor to downsample in time
        ffactor  -- Factor to downsample in frequecy (must be factor of nchans)
        filename -- File to write downsampled data to.

        Returns: Filterbank instance (from output file)
        """
        if filename is None:
            filename = "%s_f%d_t%d.fil"%(self.basename,ffactor,tfactor)
        self.f.seek(self.hdrlen)
        if not self.info["nchans"]%ffactor == 0:
            raise ValueError,"Bad frequency factor given"
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
        foldBuffer  = Buffer(nbins*nints*nbands,C.c_float)
        countBuffer = Buffer(nbins*nints*nbands,C.c_int)
        readBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readSkipBack(gulp,maxDelay)
        for nsamps,ii in passPlan.makePass():
            lib.foldFil(readBuffer.Cbuffer, foldBuffer.Cbuffer, countBuffer.Cbuffer,
                        delayPointer, maxDelay, C.c_double(self.info["tsamp"]),
                        C.c_double(period), gulp, self.info["nchans"], nbins,
                        nints, nbands, ii*(gulp-maxDelay))
        parentInfo = self.info.copy()
        return FoldedData(parentInfo,foldBuffer,countBuffer,
                          period,dm,nbins,nints,nbands)
        
        
    
    def getStatistics(self,window=10001):        

        self.f.seek(self.hdrlen)
        gulp = 2*window

        lastread = (self.info["nsamples"]-gulp)%(gulp-window)
        nreads = (self.info["nsamples"]-gulp)//(gulp-window)
        
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

        
        # Deal with first read
        self.f.cread(readBuffer,gulp*self.info["nchans"])
        getStats(readBuffer.Cbuffer, meansBuffer.Cbuffer,
                     bpassBuffer.Cbuffer, stdevBuffer.Cbuffer,
                     writeBuffer.Cbuffer,maximaBuffer.Cbuffer,
                     minimaBuffer.Cbuffer,self.info["nchans"],gulp,window,0)
        outFile.cwrite(writeBuffer)
        
        for ii in xrange(nreads):
            sys.stdout.write("%d%%\r"%(100*ii/nreads))
            sys.stdout.flush()
            self.f.seek(-window*self.info["nchans"]*C.sizeof(readBuffer.dtype),os.SEEK_CUR)
            self.f.cread(readBuffer,gulp*self.info["nchans"])

            getStats(readBuffer.Cbuffer, meansBuffer.Cbuffer,
                     bpassBuffer.Cbuffer, stdevBuffer.Cbuffer,
                     writeBuffer.Cbuffer,maximaBuffer.Cbuffer,
                     minimaBuffer.Cbuffer,self.info["nchans"],gulp,window,1)

            outFile.cwrite(writeBuffer,
                           nunits = (gulp-window)*self.info["nchans"],
                           offset = self.info["nchans"]*window)

        self.f.seek(-window*self.info["nchans"]*C.sizeof(readBuffer.dtype),os.SEEK_CUR)
        self.f.cread(readBuffer,(lastread+window)*self.info["nchans"])
        getStats(readBuffer.Cbuffer, meansBuffer.Cbuffer, 
                 bpassBuffer.Cbuffer, stdevBuffer.Cbuffer,
                 writeBuffer.Cbuffer,maximaBuffer.Cbuffer,
                 minimaBuffer.Cbuffer,self.info["nchans"],lastread,window,1)

        outFile.cwrite(writeBuffer, nunits=lastread*self.info["nchans"],
                       offset=self.info["nchans"]*window)

        stdevBuffer.Ndarray = numpy.sqrt(stdevBuffer.Ndarray/self.info["nsamples"])
        
        infoFile = open("%s_RM.info"%(self.basename),"w+")
        info = {"sigma":stdevBuffer.Ndarray,
                "bandpass":bpassBuffer.Ndarray,
                "maxima":maximaBuffer.Ndarray,
                "minima":minimaBuffer.Ndarray}
        cPickle.dump(info,infoFile)
        infoFile.close()
        return Filterbank(outFile.name),readBuffer
        
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


class ReadPlan:
    def __init__(self,ObjInst,readBuffer):
        self.I = ObjInst
        self.f = self.I.f
        self.readBuffer = readBuffer
        self.f.seek(self.I.hdrlen)
        self.nsamps = self.I.info["nsamples"]
        self.nchans = self.I.info["nchans"]
        self.nreads = 0
        self.blocks = []

    def readOnce(self,readsamps=512):
        nreads = self.nsamps//readsamps
        lastread = self.nsamps%readsamps
        self.blocks = [ (ii,readsamps*self.nchans,0) for ii in range(nreads) ]
        self.blocks.append( (nreads,lastread*self.nchans,0) )
        self.nreads = nreads+1
        self.printInfo()

    def readSkipBack(self,readsamps=512,skipback=0):
        skipback = abs(skipback)
        if skipback >= readsamps:
            raise ValueError,"readsamps must be > skipback value"
        nreads = self.nsamps//(readsamps-skipback)-1
        lastread = (self.nsamps%(readsamps-skipback))+skipback
        self.blocks = [(ii,readsamps*self.nchans,-skipback*self.nchans) for ii in range(nreads)]
        self.blocks.append((nreads,lastread*self.nchans,0))
        self.nreads = nreads+1
        self.printInfo()

    def printInfo(self):
        print
        print "Filterbank reading plan:"
        print "------------------------"
        print "Number of reads:      ",self.nreads
        print "Nsamps of first read: ",self.blocks[0][1]/self.nchans
        print "Nsamps of last read:  ",self.blocks[-1][1]/self.nchans
        print "Nsamps to skip back:  ",-1*self.blocks[0][2]/self.nchans

    def readCustom(self):
        self.blocks = []

    def makePass(self):
        for ii,block,skip in self.blocks:
            self.f.cread(self.readBuffer,block)
            self.f.seek(skip,os.SEEK_CUR)
            yield block/self.nchans,ii

              
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
        outFile.write(self.headers[0].write_header())
        clags = Buffer(self.lags.shape[0],C.c_int)
        clags.Ndarray = self.lags
        for fil,offset in zip(self.filterbanks,self.lags):
            fil.f.seek(fil.hdrlen+offset)
        minlen = min([ fil.info["nsamples"]-offset for fil,offset in zip(self.filterbanks,self.lags)])
        
        readBuffers = Buffer(gulp*self.info["nchans"],C.c_ubyte,dim=self.nfiles)
        writeBuffer = Buffer(gulp*self.info["nchans"],C.c_ubyte)
        lastread = int(minlen%gulp)
        nreads = int(minlen//gulp)
        
        for ii in xrange(nreads):
            sys.stdout.write("%d%%\r"%(100*ii/nreads))
            sys.stdout.flush()
            for jj,fil in enumerate(self.filterbanks):
                fil.f.cread(readBuffers,n=jj)
                
            lib.genMBmask(readBuffers.Cbuffer,writeBuffer.Cbuffer,
                          threshold,self.nfiles,gulp,self.info["nchans"])
            outFile.cwrite(writeBuffer)

        for ii,fil in enumerate(self.filterbanks):
            fil.f.cread(readBuffers,nunits=lastread*self.info["nchans"],n=ii)

        lib.genMBmask(readBuffers.Cbuffer,writeBuffer.Cbuffer,
                      threshold,self.nfiles,lastread,self.info["nchans"])
        outFile.cwrite(writeBuffer,lastread*self.info["nchans"])
                
from TimeSeries import TimeSeries,MultiTimeSeries
from FoldedData import FoldedData
from Bandpass import BandpassFromBuffer            
        
