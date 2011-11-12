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

void getTim(float* inbuffer,float* outbuffer,int nchans,int nsamps,int index){
  int ii,jj,val;
  for (ii=0;ii<nsamps;ii++){
    for (jj=0;jj<nchans;jj++){
      outbuffer[index+ii]+=inbuffer[(nchans*ii)+jj];
    }
  }
}

void getBpass(float* inbuffer,float* outbuffer,int nchans,int nsamps){
  int ii,jj;
  for (ii=0;ii<nsamps;ii++){
    for (jj=0;jj<nchans;jj++){
      outbuffer[jj]+=inbuffer[(nchans*ii)+jj];
    }
  }
}

void dedisperse(float* inbuffer,
		float* outbuffer,
		int* delays,
		long maxdelay,
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

void getStats(float* inbuffer,
	      float* rollingSum,
	      float* bandpass,
	      float* stdev,
	      float* outbuffer,
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
	outbuffer[(ii*nchans)+jj] = inbuffer[(ii*nchans)+jj] - (rollingSum[jj]/(ii+1));
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
      outbuffer[(ii*nchans)+jj] = inbuffer[(ii*nchans)+jj] - (rollingSum[jj]/window);
      stdev[jj] += pow(inbuffer[(ii*nchans)+jj]-(rollingSum[jj]/window),2);
      if (outbuffer[(ii*nchans)+jj]>maxbuffer[jj])
	maxbuffer[jj] = outbuffer[(ii*nchans)+jj];
      else if (outbuffer[(ii*nchans)+jj]<minbuffer[jj])
	minbuffer[jj] = outbuffer[(ii*nchans)+jj];
    }
  }
}

void to8bit(float* inbuffer,
	    unsigned char* outbuffer,
	    float* factbuffer,
	    float* plusbuffer,
	    int nsamps,
	    int nchans)
{
  int ii,jj;
  for (ii=0;ii<nsamps;ii++){
    for (jj=0;jj<nchans;jj++){
      outbuffer[(ii*nchans)+jj] = inbuffer[(ii*nchans)+jj]/factbuffer[jj] - plusbuffer[jj]; 
    }
  }
}

void flagBlock(float* inbuffer,
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



