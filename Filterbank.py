import os,numpy,pylab,sys,cPickle,time
import ctypes as C
from Header import Header,MultiHeader
from SigPyProcUtils import File,Buffer,arrayToPointer


class Filterbank(Header):
    """Class to handle filterbank files in all their glory.
    Parameters
    ----------
    \tfilename   -- name of filterbank file

    Methods
    -------
    \tcollapse     -- Collapse filterbank in frequency and/or time
    \tdedisperse   -- Simple dedispersion algorithm (fixed some sigproc bugs)
    \tdownsample   -- Downsample filterbank in frequency and/or time
    """

    def __init__(self,filename):
        """Create Filterbank instance.

        Args:
        filename -- string containing name of file to process
        """
        Header.__init__(self,filename)
        self.f = File(filename,"r",nbits=self.info["nbits"])
        infofile = "%s.info"%(self.basename)
        if os.path.isfile(infofile):
            try:
                self.stats = cPickle.load(open(infofile,"r"))
            except:
                self.stats = None
        else:
            self.stats = None
        if self.ctype == C.c_float:
            self.lib = C.CDLL("libSigPyProc32.so")
        else:
            self.lib = C.CDLL("libSigPyProc8.so")


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
        readBuffer = Buffer(blocksize*self.info["nchans"],self.ctype)
        timBuffer = Buffer(self.info["nsamples"],C.c_float)
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readOnce(blocksize)
        for nsamps,ii in passPlan.makePass():
            self.lib.getTim(readBuffer.Cbuffer,timBuffer.Cbuffer,
                        self.info["nchans"],nsamps,ii*blocksize)
        return TimeSeries(self.info.copy(),timBuffer)

    def bandpass(self,blocksize=512):
        """Collapse a filterbank in time.

        Args:
        blocksize -- number of samples to read in each gulp
        
        Return: BandpassFromBuffer instance
        """
        readBuffer = Buffer(blocksize*self.info["nchans"],self.ctype)
        bpassBuffer = Buffer(self.info["nchans"],C.c_float)
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readOnce(blocksize)
        for nsamps,ii in passPlan.makePass():
            self.lib.getBpass(readBuffer.Cbuffer,bpassBuffer.Cbuffer,
                          self.info["nchans"],blocksize)
        return BandpassFromBuffer(self.info.copy(),bpassBuffer)

    def dedisperse(self,dm,gulp=10000):
        """Dedisperse filterbank to timeseries.

        Args:
        dm        -- Dispersion measure to dedisperse to 
        gulp      -- size of block to read at a time, if chosen gulp
                     is less than maximum dispersion delay gulp is taken as 2 * max delay.

        Returns: TimeSeries instance
        """
        chanDelays = self.getDMdelays(dm)
        delayPointer = arrayToPointer(chanDelays)
        maxDelay = int(chanDelays.max())
        gulp = max(2*maxDelay,gulp)
        timLen = self.info["nsamples"]-maxDelay
        timBuffer = Buffer(timLen,C.c_float)
        readBuffer = Buffer(self.info["nchans"]*gulp,self.ctype)
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readSkipBack(gulp,maxDelay)
        for nsamps,ii in passPlan.makePass():
            self.lib.dedisperse(readBuffer.Cbuffer,timBuffer.Cbuffer,delayPointer,
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

        outFile = self.prepOutfile(filename,
                                   (("tsamp",self.info["tsamp"]*tfactor),
                                    ("nchans",self.info["nchans"]/ffactor),
                                    ("foff",self.info["foff"]*ffactor)))

        readBuffer = Buffer(tfactor*self.info["nchans"],self.ctype)
        writeBuffer = Buffer(self.info["nchans"]/ffactor,self.ctype)
        tempBuffer = Buffer(self.info["nchans"]/ffactor,C.c_int)
        self.lib.downsampleFil(self.f.f, outFile.f, readBuffer.Cbuffer,
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
        readBuffer = Buffer(self.info["nchans"]*gulp,self.ctype)
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readSkipBack(gulp,maxDelay)
        for nsamps,ii in passPlan.makePass():
            self.lib.foldFil(readBuffer.Cbuffer, foldBuffer.Cbuffer, countBuffer.Cbuffer,
                        delayPointer, maxDelay, C.c_double(self.info["tsamp"]),
                        C.c_double(period), gulp, self.info["nchans"], nbins,
                        nints, nbands, ii*(gulp-maxDelay))
        parentInfo = self.info.copy()
        return FoldedData(parentInfo,foldBuffer,countBuffer,
                          period,dm,nbins,nints,nbands)

    def dropBits32to8(self,gulp=1024,flag=0,clip=0):
        if self.stats is None:
            raise AttributeError,"First run getStatistics()"
        
        readBuffer = Buffer(self.info["nchans"]*gulp,self.ctype)
        writeBuffer = Buffer(self.info["nchans"]*gulp,C.c_ubyte)
        
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readOnce(gulp)
        
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
        chanFactor = ((chanMax-chanMin)/255.)
        chanPlus = (chanMin/chanFactor)
        outFile = self.prepOutfile("%s_8bit.fil"%(self.basename),(("nbits",8)))
        
        if flag is not 0:
            changes = (("nbits",8),("source_name","%s_mask"%(self.info["source_name"].split()[0])))
            flagFile = self.prepOutfile("%s_8bit_mask.fil"%(self.basename),changes)

        for nsamps,ii in passPlan.makePass():
            self.lib.to8bit(readBuffer.Cbuffer,writeBuffer.Cbuffer,
                            arrayToPointer(chanFactor),arrayToPointer(chanPlus),
                            nsamps,self.info["nchans"])
            outFile.cwrite(writeBuffer)
            if flag is not 0:
                lib.flagBlock(readBuffer.Cbuffer,flagBuffer.Cbuffer,
                              arrayToPointer(flagMax),arrayToPointer(flagMin),
                              nsamps,self.info["nchans"])
                flagFile.cwrite(flagBuffer)
        if flag is not 0:
            return(Filterbank(outFile.name),Filterbank(flagFile.name))
        else:
            return Filterbank(outFile.name)

    def getStatistics(self,window=10001,gulp=30003):

        if gulp < window: raise ValueError,"gulp must be > window"              

        readBuffer = Buffer(self.info["nchans"]*gulp,self.ctype)
        writeBuffer = Buffer(self.info["nchans"]*gulp,self.ctype)
        maximaBuffer = Buffer(self.info["nchans"],C.c_float)
        minimaBuffer = Buffer(self.info["nchans"],C.c_float)
        meansBuffer = Buffer(self.info["nchans"],C.c_float)
        bpassBuffer = Buffer(self.info["nchans"],C.c_float)
        stdevBuffer = Buffer(self.info["nchans"],C.c_float)

        outFile = self.prepOutfile("%s_RM.fil"%(self.basename))
        passPlan = ReadPlan(self,readBuffer)
        passPlan.readSkipBack(gulp,window)

        for nsamps,ii in passPlan.makePass():
            self.lib.getStats(readBuffer.Cbuffer, meansBuffer.Cbuffer,
                     bpassBuffer.Cbuffer, stdevBuffer.Cbuffer,
                     writeBuffer.Cbuffer,maximaBuffer.Cbuffer,
                     minimaBuffer.Cbuffer,self.info["nchans"],nsamps,window,ii)

            if ii == 0:
                outFile.cwrite(writeBuffer)
            elif ii == passPlan.nreads-1:
                outFile.cwrite(writeBuffer, nunits=nsamps*self.info["nchans"],
                               offset=self.info["nchans"]*window)
            else:
                outFile.cwrite(writeBuffer,
                               nunits = (gulp-window)*self.info["nchans"],
                               offset = self.info["nchans"]*window)

        stdevBuffer.Ndarray = numpy.sqrt(stdevBuffer.Ndarray/self.info["nsamples"])

        infoFile = open("%s_RM.info"%(self.basename),"w+")
        info = {"sigma":stdevBuffer.Ndarray,
                "bandpass":bpassBuffer.Ndarray,
                "maxima":maximaBuffer.Ndarray,
                "minima":minimaBuffer.Ndarray}
        cPickle.dump(info,infoFile)
        infoFile.close()
        return Filterbank(outFile.name)

    def prepOutfile(self,filename,headerChanges=None):
        outFile = File(filename,"w+")
        if headerChanges is not None:
            for key,val in headerChanges:
                self.modify_header(key,val)
        outFile.write(self.write_header())
        self.reset_header()
        return outFile


class ReadPlan:
    def __init__(self,ObjInst,readBuffer,multi=False):
        self.I = ObjInst
        self.multi = multi
        self.readBuffer = readBuffer
        if self.multi:
            self.nfiles = ObjInst.nfiles
            self.fils = [inst for inst in self.I.filterbanks]
            for inst in self.fils: inst.f.seek(inst.hdrlen)
        else:
            self.f = self.I.f
            self.f.seek(self.I.hdrlen)
        
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
        nreads = self.nsamps//(readsamps-skipback)
        lastread = self.nsamps-(nreads*(readsamps-skipback))
        if lastread<skipback:
            lastread += skipback
            nreads -= 1
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
        print
        
    def makePass(self):
        start = time.time()
        for ii,block,skip in self.blocks:
            sys.stdout.write("Percentage complete: %d%%\r"%(100*ii/self.nreads))
            sys.stdout.flush()
            if self.multi:
                for jj in range(self.nfiles):
                    self.fils[jj].f.cread(self.readBuffer,nunits=block,n=jj)
                    self.f.seek(skip*self.dtypesize,os.SEEK_CUR)
            else:
                self.f.cread(self.readBuffer,nunits=block)
                self.f.seek(skip*self.dtypesize,os.SEEK_CUR)
            yield block/self.nchans,ii
        print "Execution time: %f seconds     \n"%(time.time()-start)
              
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
        self.lags*=self.info["nchans"]
        return zeroDMs,self.lags

    def sumPols(self):
        outFile = self.filterbanks[0].prepOutfile("SummedPols.fil")
        readBuffers = Buffer(gulp*self.info["nchans"],self.ctype,dim=2)
        writeBuffer = Buffer(gulp*self.info["nchans"],self.ctype)
        passPlan = ReadPlan(self,readBuffers,multi=True)
        passPlan.readOnce()
        for nsamps,ii in passPlan.makePass():
            lib.sumPols(readBuffers.Cbuffer,writeBuffer.Cbuffer,nsamps,self.info["nchans"])
            outFile.cwrite(writeBuffer)
        return Filterbank(outFile)
        
    def genMBmask(self,gulp=512,threshold=4):
        outFile = self.filterbanks[0].prepOutfile("MBmask.fil")
        clags = Buffer(self.lags.shape[0],C.c_int)
        clags.Ndarray = self.lags

        for fil,offset in zip(self.filterbanks,self.lags):
            fil.f.seek(fil.hdrlen+offset)
        minlen = min([ fil.info["nsamples"]-offset for fil,offset in zip(self.filterbanks,self.lags)])
        readBuffers = Buffer(gulp*self.info["nchans"],C.c_ubyte,dim=self.nfiles)
        writeBuffer = Buffer(gulp*self.info["nchans"],C.c_ubyte)
        passPlan = ReadPlan(self,readBuffers,multi=True)
        passPlan.readOnce()
        for nsamps,ii in passPlan.makePass():
            self.lib.genMBmask(readBuffers.Cbuffer,writeBuffer.Cbuffer,
                               threshold,self.nfiles,nsamps,self.info["nchans"])
            outFile.cwrite(writeBuffer,nsamps*self.info["nchans"])

                
from TimeSeries import TimeSeries,MultiTimeSeries
from FoldedData import FoldedData
from Bandpass import BandpassFromBuffer            
        
