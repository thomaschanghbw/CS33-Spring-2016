#include "func.h"

void Func1(int c[][TSIZE], int a[][TSIZE], int b[][TSIZE])
{
	int i;
	int blockSize = 36;
	
	#pragma omp parallel for
	for(i=0;i<TSIZE;i+=blockSize)
	  {
	    int j;
	    int end;
	    if(i+blockSize>TSIZE)
	      end =TSIZE;
	    else
	      end =i+blockSize;
	    for(j=0;j<TSIZE;j+=blockSize)
	      {
		int h;
		int end2;
		if(j+blockSize>TSIZE)
		  end2 = TSIZE;
		else
		  end2 = j+blockSize;
		for(h=i; h<end; h++)
		  {
		    int l;
		    for(l=j; l<end2; l++)
		      {
			int temp = 0;
			int k;
			for(k=0; k<TSIZE; k++)
			  {
			    temp+=a[h][k]*b[k][l];
			  }
			c[h][l]+=temp;
		      }
		  }
	      }
	  }
	/*	
       	#pragma omp parallel
       	{
#pragma omp for
	  for (i=0; i<TSIZE; i++) {
	    int j;
	    //        	    #pragma omp parallel for
	    for (j=0; j<TSIZE; j++) {
	      	      int temp = 0;
	      int k;
	      //#pragma omp parallel for
	      for (k=0; k<TSIZE; k++) {
		temp+=a[i][k]*b[k][j];
	      }
	      //	      #pragma omp critical
	        c[i][j]+=temp;
	    }
      	  }
	  } */
}

void Func2(int d[][MSIZE], int c[][MSIZE])
{
  //int i;
	/*	
	#pragma omp parallel
	{
	  int id = omp_get_thread_num();
	  int numThreads = omp_get_num_threads();
	  int chunkSize = MSIZE/numThreads;
	  int end;
	  if(id == numThreads-1)
	    end = MSIZE;
	  else
	    end = id*chunkSize + chunkSize;
	  //	  #pragma omp parallel for
	  for(i=omp_get_thread_num()*chunkSize; i<end; i++)
	    {
	      int j;
	      for(j=0; j<MSIZE; ++j)
		d[i][j]=c[j][i];
	    }
	}
	*/
	  /*
	  for(i=0; i< MSIZE; i+=chunkSize)
	    {
	      int j;
	      for(j=0; j<MSIZE; j+=chunkSize)
		{
		  int k;
		  for(k=i; k<end; ++k)
		    {
		      int l;
		      for(l=j; j<end; ++j)
			{
			  d[k][l] = c[l][k];
			}
		    }
		}
	    }
	}
	  */
  /*
	int blockSize = 50;
	#pragma omp parallel
	{
	  int i;
	  blockSize = MSIZE / omp_get_num_threads();
	  int start = blockSize*omp_get_thread_num();
	  for(i=start; i<start+blockSize ; i+=blockSize)
	  {
	    int j;
	    int end;
	    if(start+2*blockSize > MSIZE)
	      end = MSIZE;
	    else
	      end = start+blockSize;
	    for(j=0; j<MSIZE; j+=blockSize)
	      {
		int end2;
		if(j+2*blockSize >MSIZE)
		   end2 = MSIZE;
		else
		  end2 = j+blockSize;
		int a;
		for(a=i; a<end; a++)
		  {
		    int b;
		    for(b=j; b<end2; b++)
		      {
			d[a][b]=c[b][a];
		      }
		  }
	      }
	  }
	} */
	/*	
       	#pragma omp parallel for
	for(i=0; i<MSIZE; i++)
	 {
	   int j;
	   //  int thread = omp_get_thread_num();
	   //int chunkSize=MSIZE/omp_get_num_threads();
	   //int end; = thread*chunkSize + chunkSize;
	   #pragma omp parallel for
	   for(j=0; j<MSIZE; j++)
	     {
	       d[i][j]=c[j][i];
	     }
	 }
*/

  int i;
  int blockSize = 24;

          #pragma omp parallel for
  for(i=0;i<MSIZE;i+=blockSize)
    {
      int j;
      int end;
      if(i+blockSize>MSIZE)
	end =MSIZE;
      else
	end =i+blockSize;
      for(j=0;j<MSIZE;j+=blockSize)
	{
	  int h;
	  int end2;
	  if(j+blockSize>MSIZE)
	    end2 = MSIZE;
	  else
	    end2 = j+blockSize;
	  for(h=i; h<end; h++)
	    {
	      int l;
	      for(l=j; l<end2; l++)
		{
		  d[h][l]=c[l][h];
		}
	    }
	}
    }
  
}

void Func3(int z[][MSIZE], int d[][MSIZE])
{
	int y;
	int near = 2;  		// The weight of neighbor
	int itself = 84; 	// The weight of center

	#pragma omp parallel for
	for(y=1; y<MSIZE-1; y++)
	  {
	    int n;

	    z[y][0] = (near*(d[y-1][0] +
			    d[y-1][0] +
			    d[y-1][1] +
			    d[y][0] +
			    d[y][1] +
			    d[y+1][0] +
			    d[y+1][0] +
			     d[y+1][1]) + itself*d[y][0])/100;
	   
	    z[y][MSIZE-1] =(near*(d[y-1][MSIZE-2] +
				 d[y-1][MSIZE-1] +
				 d[y-1][MSIZE-1] +
				 d[y][MSIZE-2] +
				 d[y][MSIZE-1] +
				 d[y+1][MSIZE-2] +
				 d[y+1][MSIZE-1] +
				 d[y+1][MSIZE-1]) +
			    itself * d[y][MSIZE-1])/100;
	    
	    for(n=1; n<MSIZE-1; n++)
	      {
		z[y][n]= (near*(d[y-1][n-1]+
			       d[y-1][n]+
			       d[y-1][n+1]+
			       d[y][n-1]+
			       d[y][n+1]+
			       d[y+1][n-1]+
			       d[y+1][n]+
				d[y+1][n+1]) + itself*d[y][n])/100;
	      }
	    
	    
	  }

	y=0;
	z[y][0] =       (near * (d[y][0] +
				 d[y][0] +
				 d[y][0] +
				 d[y][1] +
				 d[y][1] +
				 d[y+1][0] +
				 d[y+1][0] +
				 d[y+1][1]) +
			 itself * d[y][0])/100;
	
	z[y][MSIZE-1] =       (near * (d[y][MSIZE-2] +
				       d[y][MSIZE-2] +
				      d[y][MSIZE-1] +
				       d[y][MSIZE-1] +
				       d[y][MSIZE-1]+
				      d[y+1][MSIZE-2] +
				      d[y+1][MSIZE-1] +
				      d[y+1][MSIZE-1]) +
			       itself * d[y][MSIZE-1])/100;
	
	int x;
	#pragma omp parallel for
	for(x=1; x<MSIZE-1; x++)
	  {
	    z[0][x] =   (    near * (d[0][x-1] +
					    d[0][x-1] +
					    d[0][x] +
					    d[0][x+1] +
					    d[0][x+1] +
					    d[1][x-1] +
					    d[1][x] +
					    d[1][x+1]) +
			     itself * d[0][x])/100;
	  }
	
	y=MSIZE-1;
	z[y][0] =   (    near * (d[y-1][0] +
				 d[y-1][0] +
				d[y-1][1] +
				d[y][0] +
				d[y][0] +
				d[y][0] +
				d[y][1] +
      				d[y][1]) +
			 itself * d[y][0])/100;
	

	z[y][y] =       (near * (d[y-1][y-1] +
				d[y-1][y] +
				 d[y-1][y] +
				d[y][y-1] +
				d[y][y-1] +
				d[y][y] +
				d[y][y] +
			
				d[y][y]) +
			 itself * d[y][y])/100;

	#pragma omp parallel for
	for(x=1; x<MSIZE-1; x++)
	  {
	    z[y][x] =       (near * (d[y-1][x-1] +
				     d[y-1][x] +
				    d[y-1][x+1] +
				    d[y][x-1] +
				    d[y][x-1] +
				     d[y][x] +
				    d[y][x+1] +
				     d[y][x+1]) +

			     itself * d[y][x])/100;
	    
	  }
	
	/*
	  
	for (y=0; y<MSIZE; y++) {
	  int x;
	  for (x=0; x<MSIZE; x++) { */
		  /*	if (y==0) {
				if (x==0) {
					z[y][x] = 	near * d[y][x] +
						near * d[y][x+1] +
						near * d[y+1][x] +
						near * d[y+1][x+1] +
						near * d[y][x] +
						near * d[y][x+1] +
						near * d[y][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				} else if (x==MSIZE-1) {
					z[y][x] = 	near * d[y][x-1] +
						near * d[y][x] +
						near * d[y+1][x-1] +
						near * d[y+1][x] +
						near * d[y][x-1] +
						near * d[y][x] +
						near * d[y][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				} else {
				  z[y][x] = 	near * d[y][x-1] +
					  near * d[y][x+1] +
					  near * d[y+1][x-1] +
					  near * d[y+1][x+1] +
					  near * d[y][x-1] +
					  near * d[y][x+1] +
					  near * d[y][x] +
					  near * d[y+1][x] +
					  itself * d[y][x];
					  } 
					  } else */ /*if (y==MSIZE-1) {
				if (x==0) {
					z[y][x] = 	near * d[y-1][x] +
						near * d[y-1][x+1] +
						near * d[y][x] +
						near * d[y][x+1] +
						near * d[y][x] +
						near * d[y][x+1] +
						near * d[y-1][x] +
						near * d[y][x] +
						itself * d[y][x];
				} else if (x==MSIZE-1) {
					z[y][x] = 	near * d[y-1][x-1] +
						near * d[y-1][x] +
						near * d[y][x-1] +
						near * d[y][x] +
						near * d[y][x-1] +
						near * d[y][x] +
						near * d[y-1][x] +
						near * d[y][x] +
						itself * d[y][x];
				} else {
					z[y][x] = 	near * d[y-1][x-1] +
						near * d[y-1][x+1] +
						near * d[y][x-1] +
						near * d[y][x+1] +
						near * d[y][x-1] +
						near * d[y][x+1] +
						near * d[y-1][x] +
						near * d[y][x] +
						itself * d[y][x];
				}
		  } else */ /*if(y!=MSIZE-1 && y!=0){
				if (x==0) {
					z[y][x] = 	near * d[y-1][x] +
						near * d[y-1][x+1] +
						near * d[y+1][x] +
						near * d[y+1][x+1] +
						near * d[y][x] +
						near * d[y][x+1] +
						near * d[y-1][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				} else if (x==MSIZE-1) {
					z[y][x] = 	near * d[y-1][x-1] +
						near * d[y-1][x] +
						near * d[y+1][x-1] +
						near * d[y+1][x] +
						near * d[y][x-1] +
						near * d[y][x] +
						near * d[y-1][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				} else {
					z[y][x] = 	near * d[y-1][x-1] +
						near * d[y-1][x+1] +
						near * d[y+1][x-1] +
						near * d[y+1][x+1] +
						near * d[y][x-1] +
						near * d[y][x+1] +
						near * d[y-1][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				}
			}
			z[y][x]/=100;
		}
	}
*/
	
}
						
void Func4(int b[], int a[])
{
	int chuck_size=MSIZE;	 
	int array_size=VSIZE/chuck_size;
	int chuck[chuck_size];
    int i, j;

    // b[0] = a[0];
    /*
    for(i=1; i<array_size; i++)
      {
	b[i]=b[i-1]+a[i];
      }
    chuck[0]=b[array_size-1];
    #pragma omp parallel for
	for(j=1; j<chuck_size; j++) {
	  int jarray = j*array_size;
	  b[jarray]=a[jarray];
	    
	      for (i=1; i<array_size; i++) {
		b[jarray+i]=b[jarray+i-1]+a[jarray+i];
	      }
	      chuck[j]=chuck[j-1]+b[(j+1)*array_size-1];
	      for(i=0; i<array_size; i++)
		b[jarray+i]+=chuck[j-1]/(j+1);
	      
    }
    */	
	//	for(j=1; j<chuck_size; j++) {
	// chuck[j]=chuck[j-1]+chuck[j];
	//	}
	/*
	#pragma omp parallel for
	for(j=1; j<chuck_size; j++) {
	  int inc;
	  int jarray = j*array_size;
	  int chuckval = chuck[j-1]/(j+1);
	  #pragma omp parallel for
		for (inc=0; inc<array_size; inc++) {
			b[jarray+inc]+=chuckval;
		}
	}
	*/
#pragma omp parallel for 
    for(j=0; j<chuck_size; j++) {
      int k;
      int jarray=j*array_size;
      b[jarray]=a[jarray];
      for (k=1; k<array_size-1; k+=7) {
	int jarrayplusk = jarray+k;
	b[jarrayplusk]=b[jarrayplusk-1]+a[jarrayplusk];
	b[jarrayplusk+1]=b[jarrayplusk]+a[jarrayplusk+1];
	b[jarrayplusk+2]=b[jarrayplusk+1]+a[jarrayplusk+2];
	b[jarrayplusk+3]=b[jarrayplusk+2]+a[jarrayplusk+3];
	b[jarrayplusk+4]=b[jarrayplusk+3]+a[jarrayplusk+4];
	b[jarrayplusk+5]=b[jarrayplusk+4]+a[jarrayplusk+5];
	b[jarrayplusk+6]=b[jarrayplusk+5]+a[jarrayplusk+6];
	if(array_size-k <= 6) {
	  b[jarrayplusk+7]=b[jarrayplusk+6]+a[jarrayplusk+7];
	  if(array_size-k>=2)
	    {
	  b[jarrayplusk+8]=b[jarrayplusk+7]+a[jarrayplusk+8];
	  if(array_size-k>=3)
	    {
	      b[jarrayplusk+9]=b[jarrayplusk+8]+a[jarrayplusk+9];
	      if(array_size-k>=4)
		{
		b[jarrayplusk+10]=b[jarrayplusk+9]+a[jarrayplusk+10];
		if(array_size-k>=5)
		  {
		    b[jarrayplusk+11]=b[jarrayplusk+10]+a[jarrayplusk+11];
		    if(array_size-k==6)
		      b[jarrayplusk+12]=b[jarrayplusk+11]+a[jarrayplusk+12];
		  }
		}
	    }
	    }
	  break;
	}
      }
            chuck[j]=b[(j+1)*array_size-1];
    }

    for(j=1; j<chuck_size; j++) {
      chuck[j]=chuck[j-1]+chuck[j];
    }
#pragma omp parallel for
    for(j=1; j<chuck_size; j++) {
      int inc;
      int jarray=j*array_size;
      int chuckie = chuck[j-1]/(j+1);
      for (inc=0; inc<array_size; inc++) {
	b[jarray+inc]+=chuckie;
      }
    }
    /*
    chuck[0]=b[array_size-1];
    for(j=1; j<chuck_size; j++) {
      int jarray = j*array_size;
      chuck[j]=chuck[j-1]+b[jarray+array_size-1];

      int chuckieCheese=chuck[j-1]/(j+1);
      //      #pragma omp parallel for
      for(i=0; i<array_size; ++i)
	{
	  b[jarray+i]+=chuckieCheese;
	}
    } */
    /*
    for(j=1; j<chuck_size; j++) {
      for (i=0; i<VSIZE/chuck_size; i++) {
	b[j*array_size+i]+=chuck[j-1]/(j+1);
      }
      }*/
}

void Func5(int b[], int a[])
{
	int i=0, j,  stride, stride2, step;
    int temp;
	long log_N=log2(VSIZE);

	#pragma omp parallel for
	for(j=0; j<VSIZE; j+=2) {
		b[j]=a[j];
		b[j+1] = a[j] + a[j+1];
	}

	//	#pragma omp parallel for
	for(i=4; i<VSIZE; i*=2) {
	  #pragma omp parallel for
		for(j=0; j<VSIZE; j+=i) {
				b[j+i-1] = b[j+i-1] + b[j+i/2-1];
		}
	}
	
	b[VSIZE-1]=0;
	//		#pragma omp parallel for
	for(i=(log_N-1); i>=0; i--) {
	  //		stride2=(2<<i)-1;
	  //	stride=(1<<i)-1;
	  //	step=stride2+1;
	  int pow2iplus1=(int)(pow(2,i+1));
	  int pow2iminus1=(int)(pow(2, i))-1;
	  int pow2iplus1minus1=pow2iplus1-1;
#pragma omp parallel for private(temp)
	  for(j=0; j<VSIZE; j+=pow2iplus1) {
                temp=b[j+pow2iminus1];
			b[j+pow2iminus1] = b[j+pow2iplus1minus1];
			b[j+pow2iplus1minus1] = temp+b[j+pow2iplus1minus1];
		}
	}
}
