#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.141592653589793

#define W_IA 16807
#define W_IM 2147483647
#define W_AM (1.0/W_IM)
#define W_IQ 127773
#define W_IR 2836
#define W_MASK 123459876
#define W_NTAB 32
#define W_NDIV (1+(W_IM-1)/W_NTAB)
#define W_EPS 1.2e-7
#define W_RNMX (1.-W_EPS)


typedef struct {
    char name[30];
    int n[3];
    int n_fft[3];
    float ***f;
    float *data;
    float ***fft_real;
    float ***fft_im;
} FUNCTION_3;



float power(x,y)
float x,y;
{
    if(x==0 && y==0) return(1.);
    else
	if(x<=0) return(0.);
	else
	    return(exp(y*log(x)));
}

float vabs(x)
float x;
{
    if(x>0) return(x);
    else return(-x);
}

int tronq(x)
float x;
{
    int r;

    r=(int)x;
    if(r>x) r--;
    return(r);
}

float rayon(x,y)
int x,y;
{
    return(power((x*x*1.+y*y*1.),.5));
}

float rayon3(x,y,z)
int x,y,z;
{
    return(power((x*x*1.+y*y*1.+z*z*1.),.5));
}

int min2(x,y)
int x,y;
{
    if(x>y) return(y);
    else return(x);
}

int max3(x,y,z)
int x,y,z;
{
  if(x>y)
    {
      if(x>z) return(x);
      else return(z);
    }
  else
    {
      if(y>z) return(y);
      else return(z);
    }
}


void SWAP(a,b)
float *a,*b;
{
	float c;

	c=*a; *a=*b; *b=c;
}



float ran1(idum)
long *idum;
{	
	int j;
	long k;
	static long iy=0;
	static long iv[W_NTAB];
	float temp;

	if(*idum<=0 || !iy)
	{
		if(-(*idum)<1) *idum=1;
		else *idum=-(*idum);
		for(j=W_NTAB+7;j>=0;j--)
		{
			k=(*idum)/W_IQ;
			*idum=W_IA*(*idum-k*W_IQ)-W_IR*k;
			if(*idum<0) *idum+=W_IM;
			if(j<W_NTAB) iv[j]= *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/W_IQ;
	*idum=W_IA*(*idum-k*W_IQ)-W_IR*k;
	if(*idum<0) *idum+=W_IM;
	j=iy/W_NDIV;
	iy=iv[j];
	iv[j]= *idum;
	if((temp=W_AM*iy)>W_RNMX) return(W_RNMX);
	else return(temp);
}



float gausnoise(idum,sigma)
long *idum;
float *sigma;
/* note: y1 et y2 sont des realisations
independantes de la meme v.a. normale.
On perd y2 en ne retournant que y1... */
{
	float x1,x2,y1,y2,dum;

	x1=ran1(idum);
	x2=ran1(idum);

	dum=sqrt(-2.*(*sigma)*(*sigma)*log(x1));

	y1=dum*cos(2*PI*x2);
	y2=dum*sin(2*PI*x2);

	return(y1);
}


float levynoise(idum,sigma,a)
long *idum;
float *sigma;
float a;
{
    float phi,phi0,e,t,s;
    float adif,a1,a2;
    float f1,f2,f3,f4,f5,f6;

    phi=PI*(ran1(idum)-.5);
    e=-log(ran1(idum));

    if(a!=1.)
	{
	    t=a-1.+((1-a<0?-1.:1.));
	    phi0=PI*t/(2.*a);
	    adif=a*(phi-phi0);
	    a1=1./a;
	    a2=a1-1.;
	    f1=sin(adif);
	    f2=power(cos(phi),a1);
	    f3=cos(phi-adif);
	    s=(f1/f2)*power(f3/e,a2);
	    return(s);
	}
    else
	{
	    f4=.5*PI-phi;
	    f5=.5*PI*e*cos(phi);
	    f6=f5/f4;
	    s=(2./PI)*(f4*tan(phi)+log(f6));
	    return(s*exp(*sigma));
	}
}




void fourn(data,nn,ndim,isign)
float data[];
int nn[];
int ndim;
int isign;
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for(ntot=1,idim=1;idim<=ndim;idim++)
		ntot*=nn[idim];
	nprev=1;
	for(idim=ndim;idim>=1;idim--)
	{
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev<<1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for(i2=1;i2<=ip2;i2+=ip1)
		{	
			if(i2<i2rev)
				for(i1=i2;i1<=i2+ip1-2;i1+=2)
					for(i3=i1;i3<=ip3;i3+=ip2)
					{
						i3rev=i2rev+i3-i2;
						SWAP(&data[i3],&data[i3rev]);
						SWAP(&data[i3+1],&data[i3rev+1]);
					}
			ibit=ip2>>1;
			while(ibit>=ip1 && i2rev>ibit)
			{
				i2rev-=ibit;
				ibit>>=1;
			}
			i2rev+=ibit;
		}
		ifp1=ip1;
		while(ifp1<ip2)
		{
			ifp2=ifp1<<1;
			theta=isign*2*PI/(ifp2/ip1);
			wtemp=sin(.5*theta);
			wpr=-2.*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.;
			wi=0.;
			for(i3=1;i3<=ifp1;i3+=ip1)
			{
				for(i1=i3;i1<=i3+ip1-2;i1+=2)
					for(i2=i1;i2<=ip3;i2+=ifp2)
					{
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1]+=tempr;
						data[k1+1]+=tempi;
					}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev*=n;
	}
}
	
	

int  FOURIER_FORMAT(n,m)
int n,m;
{
    if(n>m/2) return(m-n);
    else return(n);
}


	

void INITIALIZE_FUNCTION_3(func,n1,n2,n3,name)
FUNCTION_3 *func;
int n1,n2,n3;
char *name;
{
	int i,j,k;

	func->n[0]=n1; func->n[1]=n2; func->n[2]=n3;
	i=1;
	while(i<=n1) i=i*2;
	func->n_fft[0]=i/2;

	i=1;
	while(i<=n2) i=i*2;
	func->n_fft[1]=i/2;

	i=1;
	while(i<=n3) i=i*2;
	func->n_fft[2]=i/2;

	func->f=(float ***)calloc(func->n[0],sizeof(float **));
	func->fft_real=(float ***)calloc(func->n_fft[0],sizeof(float **));
	func->fft_im=(float ***)calloc(func->n_fft[0],sizeof(float **));

	for(i=0;i<func->n[0];i++)
		func->f[i]=(float **)calloc(func->n[1],sizeof(float *));
	for(i=0;i<func->n[0];i++)
	    for(j=0;j<func->n[1];j++)
		func->f[i][j]=(float *)calloc(func->n[2],sizeof(float));

	for(i=0;i<func->n_fft[0];i++)
	{
		func->fft_real[i]=(float **)calloc(func->n_fft[1],sizeof(float *));
		func->fft_im[i]=(float **)calloc(func->n_fft[1],sizeof(float *));		
	}
	for(i=0;i<func->n_fft[0];i++)
	    for(j=0;j<func->n_fft[1];j++)
		{
		    func->fft_real[i][j]=(float *)calloc(func->n_fft[2],sizeof(float));
		    func->fft_im[i][j]=(float *)calloc(func->n_fft[2],sizeof(float));
		}
	func->data=(float *)calloc((func->n_fft[0]*func->n_fft[1]*func->n_fft[2]*2+1),sizeof(float));
	strcpy(func->name,name);
}


void FREE_FUNCTION_3(func)
FUNCTION_3 *func;
{
	int i,j;

	for(i=0;i<func->n[0];i++)
	  {
	    for(j=0;j<func->n[1];j++)
	      {
		free(func->f[i][j]);
		free(func->fft_real[i][j]);
		free(func->fft_im[i][j]);
	      }
	    free(func->f[i]);
	    free(func->fft_real[i]);
	    free(func->fft_im[i]);
	  }
	free(func->f);
	free(func->fft_real);
	free(func->fft_im);
	free(func->data);
	free(func->n);
	free(func->n_fft);
	free(func->name);
	free(func);
}


void PREPARE_PHYSIQUE_3(func)
FUNCTION_3 *func;
{
	int i,j,k;
	
	for(i=0;i<func->n_fft[0];i++)
		for(j=0;j<func->n_fft[1];j++)
		    for(k=0;k<func->n_fft[2];k++)
			{
			    func->data[i*func->n_fft[1]*func->n_fft[2]*2+j*func->n_fft[2]*2+k*2+1]=func->f[i][j][k];
			    func->data[i*func->n_fft[1]*func->n_fft[2]*2+j*func->n_fft[2]*2+k*2+2]=0;
			}
}


void PREPARE_SPECTRAL_3(func)
FUNCTION_3 *func;
{
	int i,j,k;

	for(i=0;i<func->n_fft[0];i++)
		for(j=0;j<func->n_fft[1];j++)
		    for(k=0;k<func->n_fft[2];k++)
			{
			    func->data[i*func->n_fft[1]*func->n_fft[2]*2+j*func->n_fft[2]*2+k*2+1]=func->fft_real[i][j][k];
			    func->data[i*func->n_fft[1]*func->n_fft[2]*2+j*func->n_fft[2]*2+k*2+2]=func->fft_im[i][j][k];
			}
}


void RECUPERE_SPECTRAL_3(func)
FUNCTION_3 *func;
{
	int i,j,k;

	for(i=0;i<func->n_fft[0];i++)
		for(j=0;j<func->n_fft[1];j++)
		    for(k=0;k<func->n_fft[2];k++)
			{
			    func->fft_real[i][j][k]=func->data[i*func->n_fft[1]*func->n_fft[2]*2+j*func->n_fft[2]*2+k*2+1];
			    func->fft_im[i][j][k]=func->data[i*func->n_fft[1]*func->n_fft[2]*2+j*func->n_fft[2]*2+k*2+2];
			}
}


void RECUPERE_PHYSIQUE_3(func)
FUNCTION_3 *func;
{
	int i,j,k;
	float mult;

	mult=(float)(func->n_fft[0]*func->n_fft[1]*func->n_fft[2]);
	for(i=0;i<func->n_fft[0];i++)
		for(j=0;j<func->n_fft[1];j++)
		    for(k=0;k<func->n_fft[2];k++)
			func->f[i][j][k]=func->data[i*func->n_fft[1]*func->n_fft[2]*2+j*func->n_fft[2]*2+k*2+1]/mult;
}

void FFT_3(func,flag)
FUNCTION_3 *func;
int flag;
{
	int nn[4];

	nn[0]=0; nn[1]=func->n_fft[0]; nn[2]=func->n_fft[1]; nn[3]=func->n_fft[2];

	if(flag==1)
	{
		PREPARE_PHYSIQUE_3(func);
		#ifdef VERBOSE
			printf("%s : %d * %d * %d values, fft made on first %d * %d * %d values\n",func->name,func->n[0],func->n[1],func->n[2],func->n_fft[0],func->n_fft[1],func->n_fft[2]);
		#endif
		fourn(func->data,nn,3,flag);
		RECUPERE_SPECTRAL_3(func);
	}
	if(flag==-1)
	{
		PREPARE_SPECTRAL_3(func);
		#ifdef VERBOSE
			printf("%s : %d * %d * %d values, inverse fft made on first %d * %d * %d values\n",func->name,func->n[0],func->n[1],func->n[2],func->n_fft[0],func->n_fft[1],func->n_fft[2]);
		#endif
		fourn(func->data,nn,3,flag);
		RECUPERE_PHYSIQUE_3(func);
	}
}





void MULTIPLY_FFT_3(f1,f2,f3)
FUNCTION_3 *f1,*f2,*f3;
{
    int i,j,k;

    for(i=0;i<f1->n_fft[0];i++)
	for(j=0;j<f1->n_fft[1];j++)
	    for(k=0;k<f1->n_fft[2];k++)
		{
		    f3->fft_real[i][j][k]=f1->fft_real[i][j][k]*f2->fft_real[i][j][k]-f1->fft_im[i][j][k]*f2->fft_im[i][j][k];
		    f3->fft_im[i][j][k]=f1->fft_real[i][j][k]*f2->fft_im[i][j][k]+f1->fft_im[i][j][k]*f2->fft_real[i][j][k];
		}
}






void NOISE_3(func,al,c1,seed)
FUNCTION_3 *func;
float al,c1;
long *seed;
{
    int i,j,k;
    float *var;

    var=(float *)calloc(1,sizeof(float));
    *var=1.;

    if(al==2.)
	for(i=0;i<func->n[0];i++)
	    for(j=0;j<func->n[1];j++)
		for(k=0;k<func->n[2];k++)
		    func->f[i][j][k]=c1*gausnoise(seed,var);
    else
	for(i=0;i<func->n[0];i++)
	    for(j=0;j<func->n[1];j++)
		for(k=0;k<func->n[2];k++)
		    func->f[i][j][k]=c1*levynoise(seed,var,al);
}




void FRACTIONNAL_FILTER_3(func,h)
/* create func the filter for h-order fractionnal integration or derivation */
FUNCTION_3 *func;
float h;
{
    int i,j,k;
    int kx,ky,kz;
    float r;
    
    for(i=0;i<func->n_fft[0];i++)
	{
	    kx=FOURIER_FORMAT(i+1,func->n_fft[0]);
	    for(j=0;j<func->n_fft[1];j++)
		{
		    ky=FOURIER_FORMAT(j+1,func->n_fft[1]);
		    for(k=0;k<func->n_fft[2];k++)
			{
			    kz=FOURIER_FORMAT(k+1,func->n_fft[2]);
			    r=rayon3(kx,ky,kz);
			    if(i!=0 || j!=0 || k!=0)
				func->fft_real[i][j][k]=power(r,h);
			}
		}
	}
}


 


void SAVE_3(func,name)
FUNCTION_3 *func;
char name[];
{
    int i,j,k;
    FILE **out;
    char **filename;

    filename=(char **)calloc(func->n[0],sizeof(char *));
    for(i=0;i<func->n[0];i++)
      filename[i]=(char *)calloc(40,sizeof(char));

    out=(FILE **)calloc(func->n[0],sizeof(FILE *));

    for(i=0;i<func->n[0];i++)
	{
	    sprintf(filename[i],"%s%d",name,i);
	    out[i]=fopen(filename[i],"w");
	    for(j=0;j<func->n[1];j++)
		{
		    for(k=0;k<func->n[2];k++)
			fprintf(out[i],"%f\t",func->f[i][j][k]);
		    fprintf(out[i],"\n");
		}
	    fclose(out[i]);
	}
}

void SAVE_3a(func,name)
FUNCTION_3 *func;
char name[];
{
    int i,j,k;
    char filename[50];
    FILE *out;

   // filename=(char **)calloc(func->n[0],sizeof(char *));
   // for(i=0;i<func->n[0];i++)
   //   filename[i]=(char *)calloc(40,sizeof(char));
  //  out=(FILE **)calloc(func->n[0],sizeof(FILE *));
   sprintf(filename,"%s",name);
   if(( out = fopen(filename,"w")) ==NULL )
   {
     printf("Error 1000\n");
     exit(0);
   }

    for(i=0;i<func->n[0];i++)
    {
    for(j=0;j<func->n[1];j++)
    {
    for(k=0;k<func->n[2];k++)
    {
    fprintf(out,"%f\n",func->f[i][j][k]);
    }
    }
    }
    fclose(out);

}


void SAVE_3b(func,name)
FUNCTION_3 *func;
char name[];
{
    int i,j,k;
    char filename[50];
    FILE *out;

   // filename=(char **)calloc(func->n[0],sizeof(char *));
   // for(i=0;i<func->n[0];i++)
   //   filename[i]=(char *)calloc(40,sizeof(char));
  //  out=(FILE **)calloc(func->n[0],sizeof(FILE *));
   sprintf(filename,"%s",name);
   if(( out = fopen(filename,"w")) ==NULL )
   {
     printf("Error 1000\n");
     exit(0);
   }

    for(k=0;k<func->n[2];k++)
    {
    for(j=0;j<func->n[1];j++)
    {
    for(i=0;i<func->n[0];i++)
    {
    fprintf(out,"%f\n",func->f[i][j][k]);
    }
    }
    }
    fclose(out);

}



main(argc,argv)
int argc;
char *argv[];
/* command line:
(1) H (Hurst exponent)
(2) Nx (length grid in x direction)
(3) Ny (length grid in y direction)
(4) Nz (length grid in z direction)
(5) seed (for random generator)
(6) name (output file name)
*/
{
	long *seed;
    	int Nx,Ny,Nz,i,j,k,n;  
    	float H;
    	char name[50],name_out[100];
    	FUNCTION_3 *noise,*filter,*fractal;

	seed=(long *)calloc(1,sizeof(long));

	sscanf(argv[1],"%f",&H);
    	sscanf(argv[2],"%d",&Nx);
    	sscanf(argv[3],"%d",&Ny);
    	sscanf(argv[4],"%d",&Nz);
    	sscanf(argv[5],"%d",seed); if(*seed>0) *seed=-*seed;
   	sscanf(argv[6],"%s",name);


    	noise=(FUNCTION_3 *)calloc(1,sizeof(FUNCTION_3)); INITIALIZE_FUNCTION_3(noise,Nx,Ny,Nz," ");
    	filter=(FUNCTION_3 *)calloc(1,sizeof(FUNCTION_3)); INITIALIZE_FUNCTION_3(filter,Nx,Ny,Nz," ");
    	fractal=(FUNCTION_3 *)calloc(1,sizeof(FUNCTION_3)); INITIALIZE_FUNCTION_3(fractal,Nx,Ny,Nz," ");
   
	NOISE_3(noise,2.,1.,seed);
	FFT_3(noise,1);
	FRACTIONNAL_FILTER_3(filter,-1.5-H);
	MULTIPLY_FFT_3(noise,filter,fractal);
	FFT_3(fractal,-1);
	//SAVE_3(fractal,name);
        SAVE_3a(fractal,name);
        //SAVE_3b(fractal,name);
}

