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

void downsampleTim(float* inbuffer,
		   float* outbuffer,
		   int factor,
		   int newLen)
{
  int ii,jj;
  for(ii=0;ii<newLen;ii++){
    for(jj=0;jj<factor;jj++)
      outbuffer[ii]+=inbuffer[(ii*factor)+jj];
    outbuffer[ii]/=(float) factor;
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


