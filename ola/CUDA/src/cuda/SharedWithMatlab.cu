#include <math.h>
#include <stdio.h>
#include <cuda_runtime.h>

#define f_(i,j) f_[(i) + (j)*(m)]
#define Z(i,j) Z[(i) + (j)*m]

__constant__ float s_H[121];	



//B is a partition of the images with dimensions thread number + patchSize
//A our extended array with padding 

__global__ void Zcalc(float const * const A,float const * const B, float *Z,float const * const H,int patchSize,float filtsigma, int m, int n)
{	
	int x = blockDim.x * blockIdx.x + threadIdx.x; 
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	int x_local=threadIdx.x;
	int y_local = threadIdx.y;
	int pad =(patchSize-1)/2;
	int dimension=blockDim.x+patchSize-1;
	extern __shared__ float Memory[];
	float *s_A=&Memory[0];//local block position in s_A
	float *g_A=&Memory[dimension*dimension];//blocks global potition in g_A
	if(x<m-2*pad && y<n-2*pad){
		s_A[x_local +y_local*dimension]=B[x_local+y_local*dimension];
		if(x_local>blockDim.x-patchSize){
			s_A[(x_local+patchSize-1) + y_local*dimension]=B[x_local+patchSize-1 +y_local*dimension];
		}
		if(y_local>blockDim.y-patchSize){
			s_A[x_local + (y_local+patchSize-1)*dimension]=B[x_local + (y_local+patchSize-1)*dimension];
		}
		if(x_local>blockDim.x-patchSize && y_local>blockDim.y-patchSize ){
			s_A[x_local+patchSize-1 + (y_local+patchSize-1)*dimension]=B[x_local+patchSize-1 + (y_local+patchSize-1)*dimension];
		}

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
		__syncthreads();
	}
	if(x<m-2*pad && y<n-2*pad){
		int counter=0;
		float temp=0,FNij=0,z_local=Z(x+pad,y+pad);
		
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
		Z[x+pad + (y+pad)*m]=z_local;
	}
}

__global__ void fCalc(float const * const A,float const * const B,float const * const Z,float const * const H,float *f_,int patchSize, float filtsigma, int m, int n){
	int x = blockDim.x * blockIdx.x + threadIdx.x; 
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	int x_local=threadIdx.x;
	int y_local = threadIdx.y;
	int pad =(patchSize-1)/2;
	int dimension=blockDim.x+patchSize-1;
	extern __shared__ float Memory[];
	float *s_A=&Memory[0];//local block position in s_A
	float *g_A=&Memory[dimension*dimension];//blocks global position g_A
	if(x<m-2*pad && y<n-2*pad){
		s_A[x_local +y_local*dimension]=B[x_local+y_local*dimension];
		if(x_local>blockDim.x-patchSize){
			s_A[(x_local+patchSize-1) + y_local*dimension]=B[x_local+patchSize-1 +y_local*dimension];
		}
		if(y_local>blockDim.y-patchSize){
			s_A[x_local + (y_local+patchSize-1)*dimension]=B[x_local + (y_local+patchSize-1)*dimension];
		}
		if(x_local>blockDim.x-patchSize && y_local>blockDim.y-patchSize ){
			s_A[x_local+patchSize-1 + (y_local+patchSize-1)*dimension]=B[x_local+patchSize-1 + (y_local+patchSize-1)*dimension];
		}
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
		__syncthreads();
	}
	if(x<m-2*pad && y<n-2*pad){
		int i,j,k,l,counter=0;
		float temp=0,FNij=0,Z_local=Z(x+pad,y+pad),f_local=f_(x+pad,y+pad);
		
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
		f_[x+pad + (y+pad)*m]=f_local;
	}
}
/*
__global__ void Zcalc(float const *const A,float const * const B,float const * const H,float *Z,int patchSize,float filtSigma,int m , int n)
{ //global coordinates
  int x=blockIdx.x*blockDim.x + threadIdx.x;
  int y=blockIdx.y*blockDim.y + threadIdx.y;
  //block coordinates
  int xblock=threadIdx.x;
  int yblock=threadIdx.y;
  extern __shared__ float Memory[];
  int dimension=blockDim.x + patchSize-1;
  
  float *s_A=&Memory[0];
  float *g_A=&Memory[dimension*dimension];
  if(x<m-patchSize+1 && y<n-patchSize+1)
  {
    s_A[xblock + yblock*dimension]=B[xblock+yblock*dimension];
    __syncthreads();
    if(xblock>blockDim.x-patchSize){
      s_A[xblock+patchSize-1+ yblock*dimension]=B[xblock+patchSize-1 + yblock*dimension];
    }
    __syncthreads();    
    if(yblock>blockDim.y-patchSize){
	  s_A[xblock+ (yblock+patchSize-1)*dimension ]=B[xblock+ (yblock+patchSize-1)*dimension ];
	}
	__syncthreads();
    if(xblock>blockDim.x-patchSize && yblock>blockDim.y-patchSize){
	    s_A[xblock+patchSize-1 + (yblock+patchSize-1)*dimension]=B[xblock+patchSize-1 + (yblock+patchSize-1)*dimension];
	}
	__syncthreads();
    //global coordinates of block
    g_A[xblock + yblock*dimension]=A[x + y*m];
    __syncthreads();
    if(xblock>blockDim.x-patchSize){
      g_A[xblock+patchSize-1 + yblock*dimension]=A[x+patchSize-1 + y*m];
    }
    __syncthreads();
    if(yblock>blockDim.y-patchSize){
	  g_A[xblock+ (yblock+patchSize-1)*dimension ]=A[x+ (y+patchSize-1)*m ];
	}
	__syncthreads();
    if(xblock>blockDim.x-patchSize && yblock>blockDim.y-patchSize){
	    g_A[xblock+patchSize-1 + (yblock+patchSize-1)*dimension]=A[x+patchSize-1 + (y+patchSize-1)*m];
	}
  }
    __syncthreads();
    patchSize=(patchSize-1)/2;

    if(x<m-2*patchSize && y<n-2*patchSize)
    {
	  int i,j,k,l,counter=0;
	  float FNij=0,temp=0;
	  float Z_local=Z(x+patchSize,y+patchSize);
	  for(i=patchSize;i<dimension-patchSize;i++){
	    for(j=patchSize;j<dimension-patchSize;j++){
		  for(k=-patchSize;k<=patchSize;k++){
		    for(l=-patchSize;l<=patchSize;l++){
		      temp=(g_A[xblock+patchSize+l + (yblock+patchSize+k)*dimension]-s_A[i+l + (j+k)*dimension])*H[counter];
		      temp=temp*temp;
		      FNij=FNij+temp;
		      counter++;
		    }
		  }
		  Z_local=Z_local+expf(- (FNij/filtSigma));
		  FNij=0;
		  counter=0;
		}
	  }
	 Z(x+patchSize,y+patchSize)=Z_local;
	}
}

__global__ void fCalc(float const *const A,float const * const B,float const * const H,float const * const Z,float *f_,int patchSize,float filtSigma,int m , int n)
{ //global coordinates
  int x=blockIdx.x*blockDim.x + threadIdx.x;
  int y=blockIdx.y*blockDim.y + threadIdx.y;
  //block coordinates
  int xblock=threadIdx.x;
  int yblock=threadIdx.y;
  extern __shared__ float Memory[];
  int dimension=blockDim.x + patchSize-1;
  
  float *s_A=&Memory[0];
  float *g_A=&Memory[dimension*dimension];
   if(x<m-patchSize+1 && y<n-patchSize+1)
  {
    s_A[xblock + yblock*dimension]=B[xblock+yblock*dimension];
    __syncthreads();
    if(xblock>blockDim.x-patchSize){
      s_A[xblock+patchSize-1+ yblock*dimension]=B[xblock+patchSize-1 + yblock*dimension];
    }
    __syncthreads();    
    if(yblock>blockDim.y-patchSize){
	  s_A[xblock+ (yblock+patchSize-1)*dimension ]=B[xblock+ (yblock+patchSize-1)*dimension ];
	}
	__syncthreads();
    if(xblock>blockDim.x-patchSize && yblock>blockDim.y-patchSize){
	    s_A[xblock+patchSize-1 + (yblock+patchSize-1)*dimension]=B[xblock+patchSize-1 + (yblock+patchSize-1)*dimension];
	}
	__syncthreads();
    //global coordinates of block
    g_A[xblock + yblock*dimension]=A[x + y*m];
    __syncthreads();
    if(xblock>blockDim.x-patchSize){
      g_A[xblock+patchSize-1 + yblock*dimension]=A[x+patchSize-1 + y*m];
    }
    __syncthreads();
    if(yblock>blockDim.y-patchSize){
	  g_A[xblock+ (yblock+patchSize-1)*dimension ]=A[x+ (y+patchSize-1)*m ];
	}
	__syncthreads();
    if(xblock>blockDim.x-patchSize && yblock>blockDim.y-patchSize){
	    g_A[xblock+patchSize-1 + (yblock+patchSize-1)*dimension]=A[x+patchSize-1 + (y+patchSize-1)*m];
	}
  }
    __syncthreads();
  patchSize=(patchSize-1)/2;
  if(x<m-2*patchSize && y<n-2*patchSize)
  {
    int i,j,k,l,counter=0;
	  float FNij=0,temp=0;
	  float Z_local=Z(x+patchSize,y+patchSize),f_local=f_(x+patchSize,x+patchSize);
	  for(i=patchSize;i<dimension-patchSize;i++){
	    for(j=patchSize;j<dimension-patchSize;j++){
		  for(k=-patchSize;k<=patchSize;k++){
		    for(l=-patchSize;l<=patchSize;l++){
		      temp=(g_A[xblock+patchSize+l + (yblock+patchSize+k)*dimension]-s_A[i+l + (j+k)*dimension])*H[counter];
		      temp=temp*temp;
		      FNij=FNij+temp;
		      counter++;
		    }
		  }
		  f_local=f_local+(1/Z_local)*(expf(- (FNij/filtSigma)))*(s_A[i+(j)*dimension]);
		  FNij=0;
		  counter=0;
		}
	  }
    f_(x+patchSize,y+patchSize)=f_local;
  }
}
*/
