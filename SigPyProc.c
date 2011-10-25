#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <time.h>
#include <Python.h>
 
/*----------------------------------------------------------------------------*/

void testarray(unsigned char* x,int lenx){

  int ii;
  for (ii=0;ii<lenx;ii++){
    printf("%d/%d\t",x[ii],ii);
  }

}

FILE* cfile(PyObject* file){
  FILE* cfile = PyFile_AsFile(file);
  return(cfile); 
}

void runningMean(float* inbuffer,
                 float* outbuffer,
                 int nsamps,
                 int window)
{
  int ii;
  double sum = 0;
  for (ii=0;ii<nsamps;ii++){
    sum += inbuffer[ii];
    if (ii<window)
      outbuffer[ii] = inbuffer[ii] - (float) sum/(ii+1);
    else {
      outbuffer[ii] = inbuffer[ii] - (float) sum/(window+1);
      sum -= inbuffer[ii-window];
    }
  }
}

void foldTim(float* buffer,
	     double* result,
	     int* counts, 
	     double tsamp,
	     double period,
	     int nsamps,
	     int nbins,
	     int nsubs)
{
  int ii,phasebin,subbint,factor;
  factor = (int) nsamps/nsubs;
  for(ii=0;ii<nsamps;ii++){
    phasebin = ((int)((ii*tsamp*nbins/period)+0.5))%nbins;
    subbint = (int) ii/factor;
    result[(subbint*nbins)+phasebin]+=buffer[ii];
    counts[(subbint*nbins)+phasebin]++;
  }
}

void ifft(float* buffer,
	  float* result,
          int size)
{
  fftwf_plan plan;
  plan = fftwf_plan_dft_c2r_1d(size, (fftwf_complex*) buffer,
			       result,FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
}

void rfft(float* buffer,
          float* result,
          int size)
{
  fftwf_plan plan;
  plan = fftwf_plan_dft_r2c_1d(size, buffer, (fftwf_complex*) result,FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
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


