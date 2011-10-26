import numpy,pylab
import ctypes as C
from Header import Header
from SigPyProcUtils import rollArray,File,Buffer
from scipy.stats.stats import chisquare
lib = C.CDLL("libSigPyProc.so")

class FoldedData(Header):
    def __init__(self,parentInfo,foldBuffer,countsBuffer,
                 period,dm,nbins,nints,nbands):
        Header.__init__(self,parentInfo)
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
                nFreqPhase[ii] = rollArray(self.freqPhase[ii],delays[ii],0)
            dmcurve[jj] = chisquare(nFreqPhase.sum(axis=0))[0]
        return dmcurve

    def optimiseP(self,loP=-0.005,hiP=0.005,Pstep=0.000005):
        periods = numpy.arange(loP,hiP,float(Pstep))
        pcurve = numpy.empty_like(periods)
        for jj,p in enumerate(periods):
            nTimePhase = numpy.empty_like(self.timePhase)
            delays = self.getPdelays(p)
            for ii in range(self.nints):
                nTimePhase[ii] = rollArray(self.timePhase[ii],delays[ii],0)
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
        
