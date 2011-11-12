#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <time.h>
#include <math.h>
#include <Python.h>
 
/*----------------------------------------------------------------------------*/

FILE* cfile(PyObject* file){
  FILE* cfile = PyFile_AsFile(file);
  return(cfile); 
}

void getTim(unsigned char* inbuffer,float* outbuffer,int nchans,int nsamps,int index){
  int ii,jj,val;
  for (ii=0;ii<nsamps;ii++){
    for (jj=0;jj<nchans;jj++){
      outbuffer[index+ii]+=inbuffer[(nchans*ii)+jj];
    }
  }
}

void getBpass(unsigned char* inbuffer,float* outbuffer,int nchans,int nsamps){
  int ii,jj;
  for (ii=0;ii<nsamps;ii++){
    for (jj=0;jj<nchans;jj++){
      outbuffer[jj]+=inbuffer[(nchans*ii)+jj];
    }
  }
}

void dedisperse(unsigned char* inbuffer,
		float* outbuffer,
		int* delays,
		int maxdelay,
		int nchans,
		int nsamps,
		int index )
{
  int ii,jj;
  
  for (ii=0;ii<(nsamps-maxdelay);ii++){
    for (jj=0;jj<nchans;jj++){
      outbuffer[index+ii] += inbuffer[(ii*nchans)+(delays[jj]*nchans)+jj];
    }
  }
}

void foldFil(unsigned char* inbuffer,
	     float* foldbuffer,
	     int* countbuffer,
	     int* delays,
	     int maxDelay,
	     double tsamp,
	     double period,
	     int nsamps,
	     int nchans,
	     int nbins,
	     int nints,
	     int nsubs,
	     int index)
{
  int ii,jj,phasebin,subband,subint,pos;
  float factor1,factor2;
  factor1 = (float) nsamps/nints;
  factor2 = (float) nchans/nsubs;

  for (ii=0;ii<(nsamps-maxDelay);ii++){
    phasebin = ((int)(((ii+index)*tsamp*nbins/period)+0.5))%nbins;
    subint = (int) ((index+ii)/factor1)%nints;
    for (jj=0;jj<nchans;jj++){
      subband = (int) jj/factor2;
      pos = (subint*nsubs*nbins)+(subband*nbins)+phasebin;
      foldbuffer[pos]+= inbuffer[(ii*nchans)+(delays[jj]*nchans)+jj];
      countbuffer[pos]++; 
    }
  }
}


void downsampleFil(FILE* infile,
		   FILE* outfile,
		   unsigned char* inbuffer,
		   unsigned char* outbuffer,
		   int* tempbuffer,
		   int tfactor,
		   int ffactor,
		   int nchans,
		   int nsamps)
{

  int ii,jj,kk,ll,iterator,val;
  int block_size = nchans*tfactor;
  int nreads = (int) (nsamps/tfactor);

  for (ii=0;ii<nreads;ii++){
    fread(inbuffer,1,block_size,infile);
    for (ll=0;ll<tfactor;ll++){
      for (kk=0;kk<nchans;kk++){
        tempbuffer[kk/ffactor] += inbuffer[(ll*nchans)+kk];
      }
    }
    for (ll=0;ll<nchans/ffactor;ll++){
      outbuffer[ll] = (unsigned char) (tempbuffer[ll]/ffactor/tfactor);
      tempbuffer[ll] = 0;
    }
    fwrite(outbuffer,1,nchans/ffactor,outfile);
  }
}

void getStats(unsigned char* inbuffer,
	      double* rollingSum,
	      double* bandpass,
	      double* stdev,
	      unsigned char* outbuffer,
	      float* maxbuffer,
	      float* minbuffer,
	      int nchans,
	      int nsamps,
	      int window,
	      int index)
{

  int ii,jj;

  if (index == 0){
    for (ii=0;ii<window;ii++){
      for (jj=0;jj<nchans;jj++){
        rollingSum[jj] += inbuffer[(ii*nchans)+jj];
        bandpass[jj] += inbuffer[(ii*nchans)+jj];
        outbuffer[(ii*nchans)+jj] = (unsigned char) (inbuffer[(ii*nchans)+jj] - (rollingSum[jj]/(ii+1)))+127;
        stdev[jj] += pow(outbuffer[(ii*nchans)+jj],2);
        if (outbuffer[(ii*nchans)+jj]>maxbuffer[jj])
          maxbuffer[jj] = outbuffer[(ii*nchans)+jj];
        else if (outbuffer[(ii*nchans)+jj]<minbuffer[jj])
          minbuffer[jj] = outbuffer[(ii*nchans)+jj];
      }
    }
  }


  for (ii=window;ii<nsamps;ii++){
    for (jj=0;jj<nchans;jj++){
      rollingSum[jj] += inbuffer[(ii*nchans)+jj];
      rollingSum[jj] -= inbuffer[((ii-window)*nchans)+jj];
      bandpass[jj] += inbuffer[(ii*nchans)+jj];
      outbuffer[(ii*nchans)+jj] = (unsigned char) (inbuffer[(ii*nchans)+jj] - (rollingSum[jj]/window))+127;
      stdev[jj] += pow(inbuffer[(ii*nchans)+jj]-(rollingSum[jj]/window),2);
      if (outbuffer[(ii*nchans)+jj]>maxbuffer[jj])
        maxbuffer[jj] = outbuffer[(ii*nchans)+jj];
      else if (outbuffer[(ii*nchans)+jj]<minbuffer[jj])
        minbuffer[jj] = outbuffer[(ii*nchans)+jj];
    }
  }
}

void flagBlock(unsigned char* inbuffer,
	       unsigned char* flagbuffer,
	       float* flagMax,
	       float* flagMin,
	       int nsamps,
	       int nchans)
{
  int ii,jj;
  for (ii=0;ii<nsamps;ii++){
    for (jj=0;jj<nchans;jj++){
      if (inbuffer[(ii*nchans)+jj] > flagMax[jj])
        flagbuffer[(ii*nchans)+jj] = 2;
      else if (inbuffer[(ii*nchans)+jj] > flagMin[jj])
        flagbuffer[(ii*nchans)+jj] = 0;
      else
        flagbuffer[(ii*nchans)+jj] = 1;
    }
  }
}


void genMBmask(unsigned char* inbuffer,
	       unsigned char* outbuffer,
	       int threshold,
	       int nfiles,
	       int nsamps,
	       int nchans)
{

  int ii,jj,kk,val,upper,lower;
  upper = nfiles+threshold;
  lower = nfiles-threshold;
    
  
  for (kk=0;kk<nfiles;kk++){
    for (ii=0;ii<nsamps*nchans;ii++){
      outbuffer[ii] += inbuffer[(kk*nsamps*nchans)+ii];
    }
  }
  for (ii=0;ii<nsamps*nchans;ii++){
    if (outbuffer[ii]>=upper || outbuffer[ii]<=lower)
      outbuffer[ii]=1;
    else
      outbuffer[ii]=0;
  }
}

