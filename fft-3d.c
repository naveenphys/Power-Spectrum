#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<fftw3.h>
#include "arrays.h"
#include "params.h"
#include "io.h"

void main()
{
  int dir; //corresponds to three velocity components
  double ****velr; //4d data array to store velocity components in r-space
  double ****velk; //4d data array to store velocity components in k-space
  int i,j,k;
  fftw_plan p;
  fftw_complex *out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*(nz/2+1));

// create a 4-d array to store velocity components
  velr = (double ****)array4d(nv,nx,ny,nz,sizeof(double));
  velk = (double ****)array4d(nv,nx,ny,(nz/2+1),sizeof(double));
  if(verbose) printarray4d(nv,nx,ny,nz,velr);

// create fftw plan
  double *in;
  in = &velr[0][0][0][0];
  p=fftw_plan_dft_r2c_3d(nx,ny,nz,in,out,FFTW_ESTIMATE);

  for(i=f1;i<f2;i++){
//   read data into the array
     read_dbl(i,velr);
     if(verbose) printarray4d(nv,nx,ny,nz,velr);
//   do fft on each component seperately
     for(dir=0;dir<nv;dir++){
        in = &velr[dir][0][0][0];
        fftw_execute(p);
	write_velk(dir,velk,out);
     }
     if(verbose) printarray4d(nv,nx,ny,nz/2+1,velk);
// create energy spectrum

// write data to output
// 0: dbl
// 1: vtk
   write_output(0,f,nv,nx,ny,nz/2+1,velk);
   write_output(1,f,nv,nx,ny,nz/2+1,velk);
  }
  freearray4d((void ****) velr);
  freearray4d((void ****) velk);
  fftw_destroy_plan(p);
  fftw_free(out);
  exit(0);
}
