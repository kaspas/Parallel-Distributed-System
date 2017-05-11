#include <math.h>
#include <stdio.h>
#include <cuda_runtime.h>


#define f(x,y) A[x + (y)* m]
#define Z(x,y) Z[x + (y)* m]
#define f_(x,y) f_[(x) + (y)*m]



__constant__ float s_H[81];	



//B is a partition of the images with dimensions thread number + patchSize
//A our extended array with padding 

__global__ void Zcalc(float const * const A, float *Z,int patchSize,float filtsigma, int m, int n)
{	
	int x = blockDim.x * blockIdx.x + threadIdx.x; 
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	int x_local=threadIdx.x;
	int y_local = threadIdx.y;
	int pad =(patchSize-1)/2;
	int dimension=blockDim.x+patchSize-1;
	extern __shared__ float Memory[];
	float z_local=0;
	float *s_A=&Memory[0];//local block position in s_A
	float *g_A=&Memory[dimension*dimension];//block's global position in g_A
	if(x<m-2*pad && y<n-2*pad){
	  z_local=Z(x+pad,y+pad);
	  g_A[x_local +y_local*dimension]=A[x+y*m];
	  if(x_local>blockDim.x-patchSize){		  
		g_A[(x_local+patchSize-1) + y_local*dimension]=A[x+patchSize-1 +y*m];
	  }
	  if(y_local>blockDim.y-patchSize){
	  g_A[x_local + (y_local+patchSize-1)*dimension]=A[x + (y+patchSize-1)*m];
	  }
	  if(x_local>blockDim.x-patchSize && y_local>blockDim.y-patchSize ){
		g_A[x_local+patchSize-1 + (y_local+patchSize-1)*dimension]=A[x+patchSize-1 + (y+patchSize-1)*m];
      }
	}
    for(int xpos=0;xpos<(m-patchSize+1);xpos=xpos+blockDim.x){
	for(int ypos=0;ypos<(n-patchSize+1);ypos=ypos+blockDim.y){	
    __syncthreads();
	if(x<m-2*pad && y<n-2*pad){
		s_A[x_local +y_local*dimension]=A[x_local+xpos+(y_local+ypos)*m];
		if(x_local>blockDim.x-patchSize){
			s_A[(x_local+patchSize-1) + y_local*dimension]=A[x_local+xpos+patchSize-1 +(y_local+ypos)*m];
		}
		if(y_local>blockDim.y-patchSize){
			s_A[x_local + (y_local+patchSize-1)*dimension]=A[x_local +xpos + (y_local+ypos+patchSize-1)*m];
		}
		if(x_local>blockDim.x-patchSize && y_local>blockDim.y-patchSize ){
			s_A[x_local+patchSize-1 + (y_local+patchSize-1)*dimension]=A[x_local+xpos+patchSize-1 + (y_local+ypos+patchSize-1)*m];
		}
	}
	__syncthreads();
	if(x<m-2*pad && y<n-2*pad){
		int counter=0;
		float temp=0,FNij=0;
		
		for(int i=pad;i<dimension-pad;i++){
			for(int j=pad;j<dimension-pad;j++){
				for(int p=-pad;p<=pad;p++){
					for(int l=-pad;l<=pad;l++){
						temp=(g_A[(x_local+pad +l)+(y_local+pad + p)*dimension]-s_A[(i+l) + (j+p)*dimension])*s_H[counter];
						FNij=FNij+temp*temp;
						counter++;
					}
				}
				z_local=z_local+expf(-(FNij/filtsigma));
				FNij=0;
				counter=0;
			}
		}
	}
}
}
  Z[x+pad + (y+pad)*m]=z_local;
}

__global__ void fCalc(float const * const A,float const * const Z,float *f_,int patchSize, float filtsigma, int m, int n){
	int x = blockDim.x * blockIdx.x + threadIdx.x; 
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	int x_local=threadIdx.x;
	int y_local = threadIdx.y;
	int pad =(patchSize-1)/2;
	int dimension=blockDim.x+patchSize-1;
	extern __shared__ float Memory[];
	float Z_local=0, f_local=0;
	float *s_A=&Memory[0];//local block position in s_A
	float *g_A=&Memory[dimension*dimension];//blocks global position g_A
    if(x<m-2*pad && y<n-2*pad){
	  Z_local=Z(x+pad,y+pad),f_local=f_(x+pad,y+pad);
      g_A[x_local +y_local*dimension]=A[x+y*m];
      if(x_local>blockDim.x-patchSize){
	    g_A[(x_local+patchSize-1) + y_local*dimension]=A[x+patchSize-1 +y*m];
	  }
	  if(y_local>blockDim.y-patchSize){
	    g_A[x_local + (y_local+patchSize-1)*dimension]=A[x + (y+patchSize-1)*m];
	  }
	  if(x_local>blockDim.x-patchSize && y_local>blockDim.y-patchSize ){
	    g_A[x_local+patchSize-1 + (y_local+patchSize-1)*dimension]=A[x+patchSize-1 + (y+patchSize-1)*m];
	  }
    }
	for(int xpos=0;xpos<(m-patchSize+1);xpos=xpos+blockDim.x){
	for(int ypos=0;ypos<(n-patchSize+1);ypos=ypos+blockDim.y){
	__syncthreads();
	if(x<m-2*pad && y<n-2*pad){
		s_A[x_local+ +y_local*dimension]=A[x_local+xpos+(y_local+ypos)*m];
		if(x_local>blockDim.x-patchSize){
			s_A[(x_local+patchSize-1) + y_local*dimension]=A[x_local+xpos+patchSize-1 +(y_local+ypos)*m];
		}
		if(y_local>blockDim.y-patchSize){
			s_A[x_local + (y_local+patchSize-1)*dimension]=A[x_local+xpos + (y_local+ypos+patchSize-1)*m];
		}
		if(x_local>blockDim.x-patchSize && y_local>blockDim.y-patchSize ){
			s_A[x_local+patchSize-1 + (y_local+patchSize-1)*dimension]=A[x_local+xpos+patchSize-1 + (y_local+ypos+patchSize-1)*m];
		}
		__syncthreads();
	}
	if(x<m-2*pad && y<n-2*pad){
		int i,j,k,l,counter=0;
		float temp=0,FNij=0;
		for(i=pad;i<dimension-pad;i++){
			for(j=pad;j<dimension-pad;j++){
				for(k=-pad;k<=pad;k++){
					for(l=-pad;l<=pad;l++){
						temp=(g_A[(x_local+pad +l)+(y_local+pad + k)*dimension]-s_A[(i+l) + (j+k)*dimension])*s_H[counter];
						FNij=FNij+temp*temp;
						counter++;
					}
				}
				f_local=f_local+(1/Z_local)*(expf(-(FNij/filtsigma)))*s_A[i+j*dimension];
				FNij=0;
				counter=0;
			}
		}
	}
}
}
  f_[x+pad + (y+pad)*m]=f_local;
}
