#include <math.h>
#include <stdio.h>
#include <cuda_runtime.h>
// Array access macros
#define f(i,j) A[(i) + (j)*(m)]
#define B(i,j) B[(i) + (j)*(m)]
#define Z(x,y) Z[(x) + (y)*(m)]
#define f_(x,y) f_[(x) + (y)*(m)]


__global__ void Zcalc(float const * const A, float *Z,float const * const H,int patchSize,float patchSigma,float fltSigma, int m, int n) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  if(x<m-(2*patchSize-1)/2 && y<n-(2*patchSize-1)/2){
    
    int i,j,k,l,counter=0;
    float FNij=0.0;
    float temp=0.0;
    patchSize=(patchSize-1)/2;
    for(i=patchSize;i<m-patchSize;i++){
	  for(j=patchSize;j<n-patchSize;j++){
		for(k=-patchSize;k<=patchSize;k++){
		  for(l=-patchSize;l<=patchSize;l++){
			temp=(f(x+patchSize+k,y+patchSize+l)-f(i+k,j+l))*H[counter];
			temp=temp*temp;
            FNij=FNij+(temp);
            counter++;
		  }
		}
	    Z(x+patchSize,y+patchSize)=Z(x+patchSize,y+patchSize)+expf(-(FNij/(fltSigma)));
	    FNij=0.0;
	    counter=0;
	  }
    }
  }
}

__global__ void fCalc(float const * const A,float const * const Z,float const * const H, float *f_,int patchSize,float patchSigma,float fltSigma, int m, int n){
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  if(x<m-(2*patchSize-1)/2 && y<n-(2*patchSize-1)/2){
	int i,j,k,l,counter=0;
    patchSize=(patchSize-1)/2;
    float FNij=0.0;
    float temp=0.0;
    float Z_local=Z(x+patchSize,y+patchSize);
    for(i=patchSize;i<m-patchSize;i++){
	  for(j=patchSize;j<n-patchSize;j++){
		for(k=-patchSize;k<=patchSize;k++){
		  for(l=-patchSize;l<=patchSize;l++){
			temp=(f(x+patchSize+k,y+patchSize+l)-f(i+k,j+l))*H[counter];
			temp=temp*temp;
            FNij=FNij+(temp);
            counter++;
		  }
		}
	    f_(x+patchSize,y+patchSize)=f_(x+patchSize,y+patchSize)+((1/Z_local)*expf(-(FNij/(fltSigma))))*f(i,j);
	    FNij=0.0;
	    counter=0;
	  }
    }
  }
}


