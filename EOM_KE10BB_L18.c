//  Gareth O'Brien 12/6/2012
// -- Include intrinsic C calls --
// mpicc -o eom EOM_KE10BB_L18.c -lm -ffast-math -lfftw3 -w
// mpirun -machinefile nodes1.pg -n 24 ./eom

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <arpa/inet.h>
#include <fftw3.h>
#include"mpi.h"
/////////////////////////////////////////////////////////////////////////////
// test parameters for real dataset

#define run_test 1

// Velocity model 
#define velocity_model 1
#define velocityStart  3000
#define velocityEnd    3000
// velocity input file model
#define input_vel  "/gpt/GEOPRO/BULK/KE/ke18_vrms8.segy"
#define NVel_ilines 2043  // velocity model length in ilines units
#define NVel_xlines 1  // velocity model number in xlines units
#define dVel_ilines 1  // velocity model iline increment
#define dVel_xlines 1  // velocity model xline increment
#define oVel_ilines 1  // velocity model iline origin
#define oVel_xlines 1  // velocity model xline origin
// image model output in (iline,xline) coordinates
#define Startiline (1)         // iline start
#define Endiline   (2043)         // iline end
#define Startxline (1)         // xline start
#define Endxline   (1)         // xline end
// image iline/xline increment
#define INCiline  (1.0)       // iline increment
#define INCxline  (1.0)       // xline increment

// Data input
#define input_data  "/gpt/MODELLING/EOM/KE/Line18/KE_10BB_2012_18_mugs_inmo.sgy"
#define output_data "/gpt/MODELLING/EOM/KE/Line18/KE_10BB_2012_18_cspVrms.sgy"

// Survey geometry
#define ox (795883.0)           // origin X
#define oy (244117.0)           // origin Y
#define angleEtoN (0.51)    // orientation of the survey angle from E to N
//#define iline_distance  (25.0)   // iline increment in metre
//#define xline_distance  (0.0)    // xline increment in metre

// CSP Gather survey
#define oxcsp (795883.0)        // origin csp X in m
#define oycsp (244117.0)        // origin csp Y in m
#define excsp (821406.0)        // end csp X in m
#define eycsp (244344.0)        // end csp Y in m
#define scaleXY  (0.001)       // scale segy positions

#define dcspX  (12.5)      // csp increment
#define dcspY  (0.0)       // csp increment
#define dh (6.25)          // csp offset increment
#define nh (1000)          // number of offsets

// geophysical constants
#define aperture (6000.0)
#define low_pass_cutoff (50.0)    // midpoint cutoff in Hz

// mpi distribution
#define image_dir "image_dir" // location of image directories
#define ntblock (20000)        // read in and process xxx traces at a time
#define number_mpi_jobsX 24   // number of image slices
#define number_mpi_jobsY 1    // number of image slices

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#define oh (0)            // number of offsets

uint32_t htonl(uint32_t hostlong);
uint16_t htons(uint16_t hostshort);
uint32_t ntohl(uint32_t netlong);
uint16_t ntohs(uint16_t netshort);

float ****alloc4d(int n1,int n2,int n3,int n4);
float ***alloc3d(int n1,int n2,int n3);
float **alloc2d(int n1,int n2);
float *alloc1d(int n1);
int *alloc1dint(int n1);
void free_array4d(float ****a,int n1,int n2,int n3,int n4);

float ntohf( const float inFloat );
void count_size_file();
void build_vel_model();
void create_lowpass_filter();
void create_sqrho_filter();
void csp(int job_id);
void write_segy_ebdic_header();
void combine_all_files();
void read_segy_header(int tracenum);

/////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265359
float sinA,cosA,csp_startX,csp_startY ;
float *anti_alias,fs,*rootrhoR,*rootrhoI;

float dtm;
int nsm,ntm;

int *xslice,*yslice;
int no_job,nx,ny;
int ns,nt,Nx,Ny,fsize,no_outtraces;
float dt,sqrho[21];
int my_rank,my_size;
int rank,size;
int dims[1],periods[1],reorder;
const  int  ndims      = 1;
MPI_Comm new_comm;
FILE *output;

/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

typedef struct {
   int tracl;
   int tracr;
   int fldr;
   int tracf;
   int ep;
   int cdp;
   int cdpt;

   short trid;
   short nvs ;
   short nhs;
   short duse;

   int offset;
   int gelev;
   int selev;
   int sdepth;
   int gdel;
   int sdel;
   int swdep;
   int gwdep;

   short scalel;
   short scalco;

   int sx;
   int sy;
   int gx;
   int gy;

   short counit;
   short wevel;
   short swevel ;
   short sut;
   short gut;
   short sstat;
   short gstat;
   short tstat ;
   short laga;
   short lagb;
   short delrt;
   short muts;
   short mute ;

   ushort ns ;
   ushort dt ;

   short gain;
   short igc;
   short igi;
   short corr;
   short sfs;
   short sfe;
   short slen;
   short styp;
   short stas;
   short stae;
   short tatype;
   short afilf;
   short afils;
   short nofilf;
   short nofils;
   short lcf;
   short hcf;
   short lcs;
   short hcs;
   short year;
   short day ;
   short hour;
   short minute;
   short sec;
   short timbas;
   short trwf;
   short grnors;
   short grnofr;
   short grnlof;
   short gaps;
   short otrav;

   int cdpx;
   int cdpy;
   int iline;
   int cline;
   int shotpt;

} headers;
/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//  MAIN CONTROL LOOP
main(int argc ,char *argv[])
{
 int job,i,j,count;
 double start,finish;
 char file_name[50];
 char command[250];

// Starting MPI calls and redefining new communicator
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

// command line arguements
  if( argc != 1 )
  {
    printf("INCORRECT USAGE - exe par_file\n");
    MPI_Finalize();  exit(1);
  }

// sorting out mpi domain 
  dims[0]    = size;
  periods[0] = 1;
  reorder    = 1;

  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &new_comm);  
  MPI_Comm_rank(new_comm, &my_rank);
  MPI_Comm_size(new_comm, &my_size);

  start=MPI_Wtime();

  if(my_rank==0)
  {
    printf("--------------------------------------------------------\n");
    printf("Equivalent Offset Migration - prototype code\n");
    printf("--------------------------------------------------------\n");
  }
//////////////////////////////////////////////////
// opening input file to determine size
  count_size_file();

 // total size of required cells for migration
  Nx=ceil( (excsp-oxcsp)/dcspX  );
  Ny=ceil( (eycsp-oycsp)/dcspY  );

  sinA=sin(angleEtoN*PI/180.0);
  cosA=cos(angleEtoN*PI/180.0);

  csp_startX =  cosA*(oxcsp-ox) + sinA*(oycsp-oy);
  csp_startY = -sinA*(oxcsp-ox) + cosA*(oycsp-oy);

  if(Nx==0)
  Nx=1;
  if(Ny==0)
  Ny=1;

  if(dcspX==0)
  Nx=1;
  if(dcspY==0)
  Ny=1;

// working out size of job volumes and no of jobs via memory requirement
  no_job=number_mpi_jobsX*number_mpi_jobsY;
  nx=(int) ceil( 1.0*Nx/number_mpi_jobsX );
  ny=(int) ceil( 1.0*Ny/number_mpi_jobsY );

  xslice=alloc1dint(no_job);
  yslice=alloc1dint(no_job);
  count=0;
  for(i=0; i<number_mpi_jobsX; i++)
  {
  for(j=0; j<number_mpi_jobsY; j++)
  {
   xslice[count]=i;
   yslice[count]=j;
   count++;
  }
  }

  if(nx==0)
  nx=1;
  if(ny==0)
  ny=1;

  if(my_rank==0)
  {
   printf("Image Size: Nx=%d Ny=%d nh=%d\n",Nx,Ny,nh);
   printf("Job Sizes:  nx=%d ny=%d nh=%d no_jobs=%d\n",nx,ny,nh,no_job);
   printf("Slices: jobsX=%d jobsY=%d\n",number_mpi_jobsX,number_mpi_jobsY);
   printf("Memory per job = %f MBytes \n",nx*ny*nh*ns*4.0/1024.0/1024.0);
   printf("Memory on input block = %f MBytes \n",ntblock*ns*4.0/1024.0/1024.0+ntblock*240.0/1024.0/1024.0);
  }

  // making output directory
  if(my_rank==0)
  {
   sprintf(command,"mkdir %s",image_dir);
   system(command);
   printf("running linux command: %s\n",command);
   printf("--------------------------------------------------------\n");
  }

 // making butter low pass filter
  create_lowpass_filter();

 // making square root rho filter
  create_sqrho_filter();

 //setting up veloicty model
 build_vel_model();

 // writing segy ebdic and binary 
  if(my_rank==0)
  write_segy_ebdic_header();

  MPI_Barrier(new_comm);

///////////////////////////////////////////////////////////////////////////////
// running mpi jobs across all cores
 if(my_rank==0)
 {
  printf("finished initialisation after %lf\n",my_rank,(MPI_Wtime()-start)/60);
  printf("--------------------------------------------------------\n");
  printf("starting building csp gathers \n");
  printf("--------------------------------------------------------\n");
 }

// main loop to generate the csp gathers
 if(run_test==1)
 {
  for(job=my_rank; job < no_job; job = job + my_size)
  {
   csp(job);
  }

  MPI_Barrier(new_comm);
  if(my_rank==0)
  combine_all_files();
 }

///////////////////////////////////////////////////////////////////////////////
// done
 if(my_rank==0)
 {
 printf("finished opening output file after %lf\n",my_rank,(MPI_Wtime()-start)/60);
 printf("--------------------------------------------------------\n");
 }

//////////////////////////////////////////////////

 finish=MPI_Wtime();
 printf("finished on node %d after %lf minutes\n",my_rank,(finish-start)/60);
 MPI_Finalize();


}          // End of main


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void count_size_file()
{
  unsigned long length;
  FILE *input;
  char file_name[50];
  ushort ns1,dt1;

  sprintf(file_name,"%s",input_data);
  if(( input = fopen(file_name,"r")) == NULL )
  {
    printf("Error 6612345: opening %s on %d\n",file_name,my_rank);
    MPI_Finalize();
    exit(1);
  }

// this may need changing for large file to determine the correct size
 fseek(input, 0, SEEK_END);
 length=ftell(input);
 fseek(input, 3200 + 400 + 114, SEEK_SET);

 fread(&ns1 , sizeof(ushort), 1, input);  ns1=htons(ns1);
 fread(&dt1 , sizeof(ushort), 1, input);  dt1=htons(dt1);

 ns= (int) ns1;
 dt= dt1/1000000.0;
 nt= ((length-3600)/(240+ns1*4));

 if(my_rank==0)
 {
  //printf("Number of samples = %lu\n" ,ns1);
  //printf("Time step in ms = %lu\n" ,dt1);
  //printf("Size of the file = %lu\n" ,length);
  //printf("Number of traces = %lu\n" , (length-3600)/(240+ns1*4));
  printf("Input File %s\n" ,file_name);
  printf("Size of the file = %lu bytes\n" ,length);
  printf("Number of samples = %d\n" ,ns);
  printf("Time step in s = %f\n" ,dt);
  printf("Number of traces = %d\n" , nt);
  printf("--------------------------------------------------------\n");
 }

 fclose(input);

 fsize=length;

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void write_segy_ebdic_header(char *argv[])
{
  short messs;
  int count,mess;
  int xc,yc,ih,xc1,yc1;
  char file_name[50];
  char ebdic[3200];
  headers trace;
  FILE *out;

  printf("CSP Rectangle array of size (%d,%d,%d)\n",Nx,Ny,nh);
  printf("--------------------------------------------------------\n");

  printf("Number of samples = %d\n" ,ns);
  printf("Time step in s = %f\n" ,dt);
  printf("Number of traces = %d\n" ,Nx*Ny*nh);
  printf("... writing output front headers\n");

///////////////////////////////////////////////////////////////////////
  /// opening edbic file
  sprintf(file_name,"%s/ebdic_header.bin",image_dir);
  if(( out = fopen(file_name,"wb")) == NULL )
  {
    printf("Error 34975765: opening %s on %d\n",file_name,my_rank);
    MPI_Finalize();
    exit(1);
  }

  // write to edbic here
  sprintf(ebdic,"outputting csp gathers from prototype code\n");
  sprintf(ebdic,"written by G. O'Brien\n");
  fwrite(ebdic, sizeof(char), 3200, out);

  fclose(out);

///////////////////////////////////////////////////////////////////////
  /// write out to binary header
  sprintf(file_name,"%s/binary_header.bin",image_dir);
  if(( out = fopen(file_name,"wb")) == NULL )
  {
    printf("Error 44975765: opening %s on %d\n",file_name,my_rank);
    MPI_Finalize();
    exit(1);
  }

  for(count=0; count<400; count=count+4)
  {
   mess=htonl(0.0);
   fwrite(&mess, sizeof(int), 1, out);
  }
  rewind(out);

  // write to binary here
  fseek(out, 12 , SEEK_SET);
  mess=htonl(Nx*Ny*nh);
  fwrite(&messs, sizeof(short), 1, out);

  fseek(out, 16 , SEEK_SET);
  messs=htons(dt*1000000);
  fwrite(&messs, sizeof(short), 1, out);

  fseek(out, 20 , SEEK_SET);
  messs=htons(ns);
  fwrite(&messs, sizeof(short), 1, out);

  fseek(out, 22 , SEEK_SET);
  messs=htons(ns);
  fwrite(&messs, sizeof(short), 1, out);

  fseek(out, 24 , SEEK_SET);       // format
  messs=htons(5);
  fwrite(&messs, sizeof(short), 1, out);

  messs=htons(1);
  fwrite(&messs, sizeof(short), 1, out);

  messs=htons(1);
  fwrite(&messs, sizeof(short), 1, out);

  messs=htons(1);
  fwrite(&messs, sizeof(short), 1, out);

  fclose(out);

  printf("--------------------------------------------------------\n");

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

void combine_all_files()
{
  int i,j,k;
  char command[250];
/*
// building output file
  printf(".. running script to combine all su files and delete them\n\n");

  sprintf(command,"cat %s/ebdic_header.bin %s/binary_header.bin > %s",image_dir,image_dir,output_data);
  printf("%s\n",command);
  system(command);

  for(i=0; i<no_job; i++)
  {
   sprintf(command,"cat %s/image%d.su >> %s",image_dir,i,output_data);
   printf("%s  ;",command);
   system(command);

   sprintf(command,"rm %s/image%d.su",image_dir,i);
   printf("  %s\n",command);
   system(command);

  }
  sprintf(command,"rm %s/ebdic_header.bin %s/binary_header.bin ",image_dir,image_dir);
  printf("%s\n",command);
  system(command);
*/
//////////////////////////////////////////////////////////////////////////////////////
  char file_name[100];
  FILE *out;

  sprintf(file_name,"%s/combine_script",image_dir);
  if(( out = fopen(file_name,"wb")) == NULL )
  {
    printf("Error 33454566: opening %s on %d\n",file_name,my_rank);
    MPI_Finalize();
    exit(1);
  }

  fprintf(out,"cat ebdic_header.bin binary_header.bin > %s \n",output_data);
  for(i=0; i<no_job; i++)
  {
   fprintf(out,"cat image%d.su >> %s\n",i,output_data);
  }
  fprintf(out,"rm image*.su\n");
  fprintf(out,"rm ebdic_header.bin binary_header.bin ");
  fclose(out);

  sprintf(command,"chmod a+x %s/combine_script",image_dir);
  system(command);
  //sprintf(command,"./%s/combine_script",image_dir);
  //system(command);

}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

float ntohf( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
float ****alloc4d(int n1,int n2,int n3,int n4)
{
  int i,j,k;
  float ****dummy;

 dummy=(float ****) calloc(n1,sizeof(float ***));
 if (dummy==NULL)
 {
    printf(" Could not allocate memory 1 4d\n");
    exit(1);
 }
 for(i=0; i<n1; i++)
 {
    dummy[i]=(float ***) calloc(n2,sizeof(float **));
    if (dummy[i]==NULL)
    {
      printf(" Could not allocate memory 2 4d \n");
      exit(1);
    }
  for(j=0; j<n2; j++)
  {
     dummy[i][j]=(float **) calloc(n3,sizeof(float *));
     if (dummy[i][j]==NULL)
     {
      printf(" Could not allocate memory 3 4d\n");
      exit(1);
     }
   for(k=0; k<n3; k++)
   {
     dummy[i][j][k]=(float *) calloc(n4,sizeof(float ));
     if (dummy[i][j][k]==NULL)
     {
       printf(" Could not allocate memory 4 4d\n");
       exit(1);
     }
   } //k
  }  //j
 }   //i 

   return dummy;
}  //end function

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void free_array4d(float ****a,int n1,int n2,int n3,int n4)
{
  int i,j,k;

 for(i=0; i<n4; i++)
 {
  for(j=0; j<n3; j++)
  {
    for(k=0; k<n2; k++)
    {
    free(a[i][j][k]);
    }
    free(a[i][j]);
  }
    free(a[i]);
 }
 free(a);

}
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

float ***alloc3d(int n1,int n2,int n3)
{
  int i,j;
  float ***dummy;

 dummy=(float ***) calloc(n1,sizeof(float **));
 if (dummy==NULL)
 {
    printf(" Could not allocate memory 1 3d\n");
    exit(1);
 }
 for(i=0; i<n1; i++)
 {
   dummy[i]=(float **) calloc(n2,sizeof(float *));
   if (dummy[i]==NULL)
   {
      printf(" Could not allocate memory 2 3d\n");
      exit(1);
   }
   for(j=0; j<n2; j++)
   {
    dummy[i][j]=(float *) calloc(n3,sizeof(float ));
    if (dummy[i][j]==NULL)
    {
      printf(" Could not allocate memory 3 3d\n");
      exit(1);
    }
   }
}

  return dummy;

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


float **alloc2d(int n1,int n2)
{
  int i;
  float **dummy;

 dummy=(float **) calloc(n1,sizeof(float *));
 if (dummy==NULL)
 {
    printf(" Could not allocate memory 1 v_2d \n");MPI_Finalize();
    exit(1);
 }
 for(i=0; i<n1; i++)
 {
   dummy[i]=(float *) calloc(n2,sizeof(float ));
   if (dummy[i]==NULL)
   {
      printf(" Could not allocate memory 2 v_2d\n");MPI_Finalize();
      exit(1);
   }
 }

 return dummy;

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
float *alloc1d(int n1)
{
 float *dummy;

 dummy=(float *) calloc(n1,sizeof(float ));
 if (dummy==NULL)
 {
    printf(" Could not allocate memory v_1d \n");MPI_Finalize();
    exit(1);
 }
  return dummy;

}
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

int *alloc1dint(int n1)
{
 int *dummy;

 dummy=(int *) calloc(n1,sizeof(int ));
 if (dummy==NULL)
 {
    printf(" Could not allocate memory i_1d \n");MPI_Finalize();
    exit(1);
 }
  return dummy;

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

void create_sqrho_filter()
{
 int i;

 fs=2.0/dt/ns;

 rootrhoR=alloc1d(ns);
 rootrhoI=alloc1d(ns);

 for(i=0;i<=(int)ns/2; i++)
 {
   rootrhoR[i]=sqrt( i*fs/2.0 );
   rootrhoI[i]=sqrt( i*fs/2.0 );
 }
 for(i=(int)ns/2;i<ns; i++)
 {
   rootrhoR[i]=rootrhoR[ns-i];
   rootrhoI[i]=rootrhoI[ns-i];
 }

 sqrho[0]=-0.0010;
 sqrho[1]=-0.0030;
 sqrho[2]=-0.0066;
 sqrho[3]=-0.0085;
 sqrho[4]=-0.0060;
 sqrho[5]=-0.0083;
 sqrho[6]=-0.0107 ;
 sqrho[7]=-0.0164;
 sqrho[8]=-0.0103;
 sqrho[9]=-0.0194;
 sqrho[10]=-0.0221 ;
 sqrho[11]=-0.0705;
 sqrho[12]=0.0395;
 sqrho[13]=-0.2161;
 sqrho[14]=-0.3831;
 sqrho[15]=0.5451;    // zero time
 sqrho[16]=0.4775;
 sqrho[17]=-0.1570;
 sqrho[18]=0.0130;
 sqrho[19]=0.0321;
 sqrho[20]=-0.0129;

}

void create_lowpass_filter()
{
 int i;
 float fmax,mm;
 char file_name[50];
 FILE *out;

 anti_alias=alloc1d(ns);
 fs=2.0/dt/ns;
 fmax=low_pass_cutoff;


 for(i=0;i<=(int)ns/2; i++)
 {
   anti_alias[i]=1.0;
 }

 for(i=0;i<=(int)ns/2; i++)
 {
   anti_alias[i]=1.0/sqrt( 1.0 + (i*fs)*(i*fs)*(i*fs)*(i*fs)/fmax/fmax/fmax/fmax );
 }
 for(i=(int)ns/2;i<ns; i++)
 {
   anti_alias[i]=anti_alias[ns-i];
 }

  if(my_rank==0)
  {
     sprintf(file_name,"%s/butterfilter.txt",image_dir);
     if(( out = fopen(file_name,"w")) == NULL )
     {
       printf("Error 9742945: opening input file %s on node %d\n",file_name,my_rank);
       MPI_Finalize();
       exit(1);
     }
     for(i=0;i<ns; i++)
     {
       fprintf(out,"%f %f\n",i*fs,anti_alias[i]);
     }
     fclose(out);
  }

}

float ***vel_model;

void build_vel_model()
{
 unsigned long length;
 ushort ns1,dt1;
 int t,it,it1,i,j,iline,xline;
 char ebdic[3600];
 char file_name[50];
 FILE *input,*inputi;
 FILE *out;


if(velocity_model==0)
{
    vel_model = alloc3d(1,1,ns);

    for(i=0;i<ns; i++)
    {
     vel_model[0][0][i] = velocityStart +  (velocityEnd-velocityStart)*(i/ns);
    }
    for(i=0;i<ns; i++)
    vel_model[0][0][i]=vel_model[0][0][i]*vel_model[0][0][i];
}

if(velocity_model==1)
{
///////////////////////////////////////////////////////////////////////////////////////
// open input file
  sprintf(file_name,"%s",input_vel);
  if(( input = fopen(file_name,"rb")) == NULL )
  {
    printf("Error 86542945: opening input file %s on node %d\n",file_name,my_rank);
    MPI_Finalize();
    exit(1);
  }

  if(( inputi = fopen(file_name,"rb")) == NULL )
  {
    printf("Error 87752945: opening input file %s on node %d\n",file_name,my_rank);
    MPI_Finalize();
    exit(1);
  }


// reading in sample information
 fseek(input, 0, SEEK_END);
 length=ftell(input);
 fseek(input, 3600 + 114, SEEK_SET);

 fread(&ns1 , sizeof(ushort), 1, input);  ns1=htons(ns1);
 fread(&dt1 , sizeof(ushort), 1, input);  dt1=htons(dt1);
 rewind(input);  // reset file scan to begining

 nsm= (int) ns1;
 dtm= dt1/1000000.0;
 ntm= ((length-3600)/(240+ns1*4));

 if(my_rank==0)
 {
  printf("--------------------------------------------------------\n");
  printf("Input Velocity Model %s\n" ,file_name);
  printf("Size of the file = %lu bytes\n" ,length);
  printf("Number of samples = %d\n" ,nsm);
  printf("Time step in s = %f\n" ,dtm);
  printf("Number of traces = %d\n" , ntm);
  printf("--------------------------------------------------------\n");
 }

// velocity model in (iline xline) units
  vel_model = alloc3d(NVel_ilines,NVel_xlines,nsm);

  int ntblock_vel=5000;             // read in xxx traces in one go
  int readin=ntblock_vel*(60+nsm);
  float *block_trace,*ibm;
  int *block_header;
  int len=readin;

  block_trace=alloc1d(readin);
  ibm=alloc1d(readin);
  block_header=alloc1dint(readin);

  // read in ebdic to force file to right spot
  fread(ebdic, 1, 3600, inputi);
  fread(ebdic, 1, 3600, input);

  int count=0;
  for(it=0; it<ntm+1; it=it+ntblock_vel ) // loop over blocks of data
  {
      fread(block_trace, 4,readin,input );   // read in data
      fread(block_header,4,readin,inputi);   // read in headers

    // converting ibm format to ieee for migration
    for (i=0; i<len; i++)
    {
        block_trace[i]=ntohf(block_trace[i]);
        ibm[i] = (((( (*(int*)&block_trace[i]) >> 31) & 0x01) * (-2)) + 1) *
                 ((float)((*(int*)&block_trace[i]) & 0x00ffffff) / 0x1000000) *
                 ((float)pow(16, (( (*(int*)&block_trace[i]) >> 24) & 0x7f) - 64));
    }

    for(it1=0; it1<ntblock_vel; it1++ ) // loop over indivdual traces in blocks
    {
      if(it1+ntblock_vel*count<ntm)
      {
       // read in iline and xline
       //iline = htonl(block_header[47 + (60+nsm)*it1]);
       //xline = htonl(block_header[48 + (60+nsm)*it1]);

   // switched xline and iline due to 2D line and issue with segy file read in iline and xline
       xline = htonl(block_header[47 + (60+nsm)*it1]);
       iline = htonl(block_header[48 + (60+nsm)*it1]);

       i=(iline-oVel_ilines)/dVel_ilines;   j=(xline-oVel_xlines)/dVel_xlines;

       if(i>=0 && i<NVel_ilines && j>=0 && j<NVel_xlines)
       {
          for(t=0; t<nsm; t++) // time loop
          {
           vel_model[i][j][t] = ibm[60 + t + (60+nsm)*it1] ; // reading in rms velocity
          }
       }

//if(my_rank==0 && iline==300 && xline==1)
//printf(" %d %d | %d %d | %f \t%d %d %d\n",iline,xline,i,j,vel_model[i][j][107],it1,it,ntm);
        }
       } // ntblock loop
    count++;
  }  // trace loop


// write out velocity model to vtk file to QC
  int mx,mi,mt;
 /* if(my_rank==0)
  {
   sprintf(file_name,"%s/migration_velocity.vtk",image_dir);
   if(( out = fopen(file_name,"w")) == NULL )
   {
    printf("Error 66642945: opening input file %s on node %d\n",file_name,my_rank);
    MPI_Finalize();
    exit(1);
   }

    i=0;mx=0; mi=0; mt=0;
    for(t=0; t<nsm; t=t+2)
    {
     mt=mt+1;mx=0;
     for(xline=0; xline<NVel_xlines; xline=xline+4)
     { 
      mx=mx+1; mi=0;
      for(iline=0; iline<NVel_ilines; iline=iline+4)
      {
       i=i+1; mi=mi+1;
      }
     }
    }

    fprintf(out,"# vtk DataFile Version 3.0\n");
    fprintf(out,"Velocity output\n");
    fprintf(out,"ASCII\n"); //fprintf(out,"BINARY\n");
    fprintf(out,"DATASET STRUCTURED_POINTS\n");
    fprintf(out,"DIMENSIONS %d %d %d\n",mi,mx,mt);
    fprintf(out,"ASPECT_RATIO 1 1 1\n");
    fprintf(out,"ORIGIN %d %d 0\n",oVel_ilines,oVel_xlines);
    fprintf(out,"POINT_DATA %d\n",mi*mx*mt);
    fprintf(out,"SCALARS Velocity float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(t=0; t<nsm; t=t+2)
    {
    for(xline=0; xline<NVel_xlines; xline=xline+4)
    {
    for(iline=0; iline<NVel_ilines; iline=iline+4)
    {
     fprintf(out,"%f \n",vel_model[iline][xline][t]);
    }
    }
    }
   fclose(out);
  }*/

 /* if(my_rank==0)
  {
   xline=0;  iline=1000;
   for(t=0; t<nsm; t++)
   printf(" %d  %d  %f  %f\n",iline,xline,t*dt,vel_model[iline][xline][t]);
  }*/


  for(iline=0; iline<NVel_ilines; iline++)
  {
  for(xline=0; xline<NVel_xlines; xline++)
  {
  for(t=0; t<nsm; t++)
  {
    if(vel_model[iline][xline][t]==0)
    printf("ERROR 11111111 - %d,%d,%d\n",iline,xline,t);

    vel_model[iline][xline][t]=vel_model[iline][xline][t]*vel_model[iline][xline][t];
  }
  }
  }

  fclose(input);
  fclose(inputi);
 }

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


void csp(int job_id)
{
 short messs;
 int it,it1,header[60];
 int xc,yc,ih,itf;
 int t,off,csp,mess,xc1,yc1;
 float scale1,scale2;
 float sxt,syt,gxt,gyt,scaleT;
 float data[ns],cspX,cspY,dist;
 float x2,h2,cmpX,cmpY,crossterm2,eq_offset;
 char ebdic[3600];
 char file_name[50];
 headers trace;
 FILE *input1,*inputi;
 FILE *out;

//////////////////////////////////////////////////
// open input file
  sprintf(file_name,"%s",input_data);
  if(( input1 = fopen(file_name,"rb")) == NULL )
  {
    printf("Error 12345: opening input file %s on node %d\n",input_data,my_rank);
    MPI_Finalize();
    exit(1);
  }

// open input file
  sprintf(file_name,"%s",input_data);
  if(( inputi = fopen(file_name,"rb")) == NULL )
  {
    printf("Error 11274345: opening input file %s on node %d\n",input_data,my_rank);
    MPI_Finalize();
    exit(1);
  }
  printf("...processing job %d of %d on node %d nt=%d \n",job_id,no_job,my_rank,nt);

//////////////////////////////////////////////////
// zeroing all trace headers
 trace.tracr=htonl(0);
 trace.fldr=htonl(0);
 trace.tracf=htonl(0);
 trace.ep=htonl(0);
 trace.cdp=htonl(0);
 trace.cdpt=htonl(0);

 trace.trid=htons(0);
 trace.nvs=htons(0);
 trace.nhs=htons(0);
 trace.duse=htons(0);

 trace.offset=htonl(0);
 trace.gelev=htonl(0);
 trace.selev=htonl(0);
 trace.sdepth=htonl(0);
 trace.gdel=htonl(0);
 trace.sdel=htonl(0);
 trace.swdep=htonl(0);
 trace.gwdep=htonl(0);

 trace.scalel=htons(0);
 trace.scalco=htons(0);

 trace.sx=htonl(13);
 trace.sy=htonl(14);
 trace.gx=htonl(15);
 trace.gy=htonl(16);

 trace.counit=htons(0);
 trace.wevel=htons(0);
 trace.swevel=htons(0);
 trace.sut=htons(0);
 trace.gut=htons(0);
 trace.sstat=htons(0);
 trace.gstat=htons(0);
 trace.tstat=htons(0);
 trace.laga=htons(0);
 trace.lagb=htons(0);
 trace.delrt=htons(0);
 trace.muts=htons(0);
 trace.mute=htons(0);

 trace.ns=htons(ns);
 trace.dt=htons(dt*1000000);

 trace.gain=htons(0);
 trace.igc=htons(0);
 trace.igi=htons(0);
 trace.corr=htons(0);
 trace.sfs=htons(0);
 trace.sfe=htons(0);
 trace.slen=htons(0);
 trace.styp=htons(0);
 trace.stas=htons(0);
 trace.stae=htons(0);
 trace.tatype=htons(0);
 trace.afilf=htons(0);
 trace.afils=htons(0);
 trace.nofilf=htons(0);
 trace.nofils=htons(0);
 trace.lcf=htons(0);
 trace.hcf=htons(0);
 trace.lcs=htons(0);
 trace.hcs=htons(0);
 trace.year=htons(0);
 trace.day =htons(0);
 trace.hour=htons(0);
 trace.minute=htons(0);
 trace.sec=htons(0);
 trace.timbas=htons(0);
 trace.trwf=htons(0);
 trace.grnors=htons(0);
 trace.grnofr=htons(0);
 trace.grnlof=htons(0);
 trace.gaps=htons(0);
 trace.otrav=htons(0);
 //181 - 240     Unassigned for optional info
 trace.cdpx=htonl(0);
 trace.cdpy=htonl(0);
 trace.iline=htonl(0);
 trace.cline=htonl(0);
 trace.shotpt=htonl(0);
 mess=htonl(0);
 messs=htons(0);
////////////////////////////////////////////////////////////////////////////
  float sqrt1,sqrt2,dsr_eq,dist_IR2,dist_IS2;
  float aperture2=aperture*aperture,io,jo,scale3,tc;
  float ****csp_data,****fold,t2d[ns],t2[ns];
  //float **block_trace;
  float *block_trace,slopeHe;
  int *block_header,i,j,talpha;
  int readin=(ntblock)*(60+ns);
  int T,count=0;

// awful coding reading in to memory twice from same file for different formats
  csp_data=alloc4d(nx,ny,nh,ns);
  //block_trace=alloc2d(ntblock,60+ns);
  block_trace=alloc1d(readin);
  block_header=alloc1dint(readin);
  fold=alloc4d(nx,ny,nh,ns);
 
  for(t=0; t<ns; t++)
  {
    t2d[t]=sqrt(t*dt)*t*dt;
    t2[t]=0.25*t*dt*t*dt;
  }

  for(xc=0; xc<nx; xc++)
  {
  for(yc=0; yc<ny; yc++)
  {
  for(ih=0; ih<nh; ih++)
  {
   for(t=0; t<ns; t++)
   {
     csp_data[xc][yc][ih][t]=0.0;   fold[xc][yc][ih][t]=1.0;
   }
  }
  }
  }
  // read in ebdic to force file to right spot
  fread(ebdic, 1, 3600, input1);
  fread(ebdic, 1, 3600, inputi);

  for(it=0; it<nt+1; it=it+ntblock ) // nt+1 loop over blocks of data
  {

    fread(block_trace,4,readin,input1);   // read in data
    fread(block_header,4,readin,inputi);  // read in headers

    // converting ibm format to ieee for migration
    for (i=0; i<readin; i++)
    {
        block_trace[i] = ntohf(block_trace[i]);
        block_trace[i] = ( (((( (*(int*)&block_trace[i]) >> 31) & 0x01) * (-2)) + 1) *
                 ((float)((*(int*)&block_trace[i]) & 0x00ffffff) / 0x1000000) *
                 ((float)pow(16, (( (*(int*)&block_trace[i]) >> 24) & 0x7f) - 64)) );
    }

    for(it1=0; it1<ntblock; it1++ ) // loop over indivdual traces in blocks
    {
      if(it1+ntblock*count<nt)
      {
       // rotate source coordinates into computational grid
       io = scaleXY*htonl(block_header[18 + (60+ns)*it1]) - ox;      jo = scaleXY*htonl(block_header[19 + (60+ns)*it1]) - oy;
       sxt =  cosA*io + sinA*jo;
       syt = -sinA*io + cosA*jo;

       // rotate receiver coordinates into computational grid
       io = scaleXY*htonl(block_header[20 + (60+ns)*it1]) - ox;      jo = scaleXY*htonl(block_header[21 + (60+ns)*it1]) - oy;
       gxt =  cosA*io + sinA*jo;
       gyt = -sinA*io + cosA*jo;

       // calculate offset and CDP coordinates in the computational grid
       h2 = 0.25*( (sxt-gxt)*(sxt-gxt) + (syt-gyt)*(syt-gyt) ) ;
       cmpX = 0.5*(sxt+gxt);
       cmpY = 0.5*(syt+gyt);

       // EOM loop
       for(xc=0; xc<nx; xc++)
       {
       for(yc=0; yc<ny; yc++)
       {
         // coordinates of image point in meters on computational grid
         cspX =  csp_startX + (xc + xslice[job_id]*nx)*dcspX;
         cspY =  csp_startY + (yc + yslice[job_id]*ny)*dcspY;
         x2 = ( (cmpX-cspX)*(cmpX-cspX) + (cmpY-cspY)*(cmpY-cspY) );
         dist=(x2 + h2);

         // distances for migration
         //dist_IS2 = (cspX-sxt)*(cspX-sxt) + (cspY-syt)*(cspY-syt);
         //dist_IR2 = (cspX-gxt)*(cspX-gxt) + (cspY-gyt)*(cspY-gyt);

         if( dist<aperture2)
         {
          // velocity model coordinates
           xc1=round( (Startiline-oVel_ilines)/dVel_ilines + (xc + xslice[job_id]*nx)*(INCiline/dVel_ilines) );
           yc1=round( (Startxline-oVel_xlines)/dVel_xlines + (yc + yslice[job_id]*ny)*(INCxline/dVel_xlines) );

           // error check on vel model
            if(xc1>NVel_ilines-1)
            xc1=NVel_ilines-1;

            if(yc1>NVel_xlines-1)
            yc1=NVel_xlines-1;

            if(velocity_model==0)
            { xc1=0;  yc1=0; }

         // talpha = round( 1.0*sqrt(x2/vel_model[xc1][yc1][0])/dt );
          talpha = round( 2.0*sqrt(x2/vel_model[xc1][yc1][0])/dt );

          for(t=talpha; t<ns; t++) // time for loop
          { 
            // double square root equation change this for vti and tti
            //sqrt1 = t2[t] + dist_IS2/vel_t[t];
            //sqrt2 = t2[t] + dist_IR2/vel_t[t];
            //dsr_eq = sqrt( sqrt1 ) + sqrt( sqrt2  ) ;
            //T= (int) floor( dsr_eq/dt );
 
            crossterm2 = (4.0*x2*h2)/(t*dt*t*dt*vel_model[xc1][yc1][t]);
            eq_offset=2.0*sqrt( x2 + h2 - crossterm2 );
            off=(int) floor((eq_offset-oh)/dh) ;

            scale2 = (eq_offset-oh)/dh - off;
            scale1 = 1.0 - scale2;

            //slopeHe=t*dt-2.0*eq_offset/sqrt(vel_t[t]);
            //tc=2.25*eq_offset/sqrt(vel_t[t]);
            //tc=2.25*eq_offset/3000.0;
            //scale3 = 0.5 + atan(200.0*(t*dt - tc)/0.080)/PI; 
            //scaleT = 1.0 - 4.0*sqrt(x2)/eq_offset;

            scaleT = 1.0;//*scale3 ;

      	    if(off<nh-1 && off>=0) //T<ns && slopeHe>0 &&
            {
              csp_data[xc][yc][off][t] = csp_data[xc][yc][off][t]     + scale1*scaleT*( block_trace[60 + t + (60+ns)*it1])/t2d[t];
                  fold[xc][yc][off][t] =     fold[xc][yc][off][t]     + 1.0 ;

              csp_data[xc][yc][off+1][t] = csp_data[xc][yc][off+1][t] + scale2*scaleT*( block_trace[60 + t + (60+ns)*it1])/t2d[t];
                  fold[xc][yc][off+1][t] =     fold[xc][yc][off+1][t] + 1.0 ;
            }
           }  // time loop
          } // aperture if

       }  // yc loop
       }  // xc loop

     } // total blocks
    }  // ntblock loop 
   count++;
  }  // trace loop

/////////////////////////////////////////////////////////////////////////////
// opening  file
  printf("...writing to file job %d on node %d it=%d\n",job_id,my_rank,it);

  sprintf(file_name,"%s/image%d.su",image_dir,job_id);
  if(( out = fopen(file_name,"wb")) == NULL )
  {
    printf("Error 87315765: opening %s on %d\n",file_name,my_rank);
    MPI_Finalize();
    exit(1);
  }

/////////////////////////////////////////////////////////////////////////////
// FFTW routines and variables
  int SIZE=ns;
  float a,b,c,d;
  float dataTr[ns];

  fftw_complex    *dataT, *fft_result, *ifft_result;
  fftw_plan       plan_forward, plan_backward;

  dataT       = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) *SIZE );
  fft_result  = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) *SIZE );
  ifft_result = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) *SIZE );
  plan_forward  = fftw_plan_dft_1d( SIZE, dataT, fft_result, FFTW_FORWARD, FFTW_ESTIMATE );
  plan_backward = fftw_plan_dft_1d( SIZE, fft_result, ifft_result, FFTW_BACKWARD, FFTW_ESTIMATE );
/////////////////////////////////////////////////////////////////////////////

  csp= nx*ny*job_id + 1;
  it = nx*ny*nh*job_id + 1;
// wont work on large files as wont right to right location
  for(xc=0; xc<nx; xc++)
  {
  for(yc=0; yc<ny; yc++)
  {
  for(ih=0; ih<nh; ih++)
  {
/*
   data[0]=0.0;
   for(t=1; t<ns; t++)
   {
     data[t] = ntohf(csp_data[xc][yc][ih][t]);
   }
*/
   //////////////////////////////////////////////
   // time domain rho filter
 /*  for(t=0; t<ns; t++)
   {
     dataTr[t]=(csp_data[xc][yc][ih][t]/fold[xc][yc][ih][t]*dt*t);
     //data[t]=ntohf(fold[xc][yc][ih][t]);
   }
   for(t=0; t<15; t++)
   data[t] =  ntohf(0.0);
   for(t=15; t<ns-5; t++)
   {
     data[t]=ntohf(sqrho[0]*dataTr[t-15] + sqrho[1]*dataTr[t-14] + sqrho[2]*dataTr[t-13] + 
                   sqrho[3]*dataTr[t-12] + sqrho[4]*dataTr[t-11] + sqrho[5]*dataTr[t-10] +
                   sqrho[6]*dataTr[t-9]  + sqrho[7]*dataTr[t-8]  + sqrho[8]*dataTr[t-7]  + 
                   sqrho[9]*dataTr[t-6]  + sqrho[10]*dataTr[t-5] + sqrho[11]*dataTr[t-4] +
                   sqrho[12]*dataTr[t-3] + sqrho[13]*dataTr[t-2] + sqrho[14]*dataTr[t-1] + 
                   sqrho[15]*dataTr[t]   + sqrho[16]*dataTr[t+1] + sqrho[17]*dataTr[t+2] + 
                   sqrho[18]*dataTr[t+3] + sqrho[19]*dataTr[t+4] + sqrho[20]*dataTr[t+5] );
   }
   for(t=ns-5; t<ns; t++)
   data[t] =  ntohf(0.0);
*/
   //////////////////////////////////////////////
   // frequency domain operations on fft_result
   t=0;  dataT[t][0] = 0.0; dataT[t][1] = 0.0;
   for(t=1; t<ns; t++)
   {
     dataT[t][0] = (csp_data[xc][yc][ih][t]/fold[xc][yc][ih][t])*dt*t;
     dataT[t][1] = 0.0;
   }
   fftw_execute( plan_forward );
   // low pass filter 
   for(t=0; t<ns; t++)
   {
     a = anti_alias[t]; b=0.0; c=fft_result[t][0]; d=fft_result[t][1];
     //a = 1.0; b=0.0; c=fft_result[t][0]; d=fft_result[t][1];
     fft_result[t][0] = a*c - b*d ;
     fft_result[t][1] = a*d + b*c ;
   }
   // sqrt rho filter
   for(t=0; t<ns; t++)
   {
     a = rootrhoR[t]; b=rootrhoI[t]; c=fft_result[t][0]; d=fft_result[t][1];
     fft_result[t][0] = a*c - b*d ;
     fft_result[t][1] = a*d + b*c ;
   }
   fftw_execute( plan_backward );
   data[0]=0.0;
   for(t=1; t<ns; t++)
   {
     data[t] = ntohf( ifft_result[t][0]/SIZE );
   }
   //////////////////////////////////////////////
   // determining trace headers for su file

   if(it<Nx*Ny*nh+1)
   {
     // determining trace headers for su file
     trace.offset=(int) ceil(oh+ih*dh);
     trace.offset=htonl(trace.offset);

     cspX =  csp_startX + (xc + xslice[job_id]*nx)*dcspX;
     cspY =  csp_startY + (yc + yslice[job_id]*ny)*dcspY;
     xc1 = cosA*cspX - sinA*cspY ;
     yc1 = sinA*cspX + cosA*cspY ;

     trace.cdpx = htonl(ox + xc1);
     trace.cdpy = htonl(oy + yc1);

     trace.tracl=htonl(it);
     trace.cdp=htonl( csp ) ;

     trace.iline=htonl(csp);
     trace.cline=htonl(1.0);

     fwrite(&trace.tracl, 4, 1, out);
     fwrite(&trace.tracr, 4, 1, out);
     fwrite(&trace.fldr,  4, 1, out);
     fwrite(&trace.tracf, 4, 1, out);
     fwrite(&trace.ep,    4, 1, out);
     fwrite(&trace.cdp,   4, 1, out);
     fwrite(&trace.cdpt,  4, 1, out);
     fwrite(&trace.trid, 2, 1, out);
     fwrite(&trace.nvs,  2, 1, out);
     fwrite(&trace.nhs,  2, 1, out);
     fwrite(&trace.duse, 2, 1, out);
     fwrite(&trace.offset, 4, 1, out);
     fwrite(&trace.gelev,  4, 1, out);
     fwrite(&trace.selev,  4, 1, out);
     fwrite(&trace.sdepth, 4, 1, out);
     fwrite(&trace.gdel,  4, 1, out);
     fwrite(&trace.sdel,  4, 1, out);
     fwrite(&trace.swdep, 4, 1, out);
     fwrite(&trace.gwdep, 4, 1, out);
     fwrite(&trace.scalel, 2, 1, out);
     fwrite(&trace.scalco, 2, 1, out);
     fwrite(&trace.sx, 4, 1, out);  
     fwrite(&trace.sy, 4, 1, out);
     fwrite(&trace.gx, 4, 1, out);  
     fwrite(&trace.gy, 4, 1, out);  
     fwrite(&trace.counit, 2, 1, out); 
     fwrite(&trace.wevel,  2, 1, out);  
     fwrite(&trace.swevel, 2, 1, out);  
     fwrite(&trace.sut,    2, 1, out); 
     fwrite(&trace.gut,    2, 1, out); 
     fwrite(&trace.sstat,  2, 1, out); 
     fwrite(&trace.gstat,  2, 1, out);  
     fwrite(&trace.tstat,  2, 1, out); 
     fwrite(&trace.laga,   2, 1, out);  
     fwrite(&trace.lagb,   2, 1, out); 
     fwrite(&trace.delrt,  2, 1, out); 
     fwrite(&trace.muts,   2, 1, out);  
     fwrite(&trace.mute ,  2, 1, out);  
     fwrite(&trace.ns , sizeof(ushort), 1, out); 
     fwrite(&trace.dt , sizeof(ushort), 1, out);  
     fwrite(&trace.gain, 2, 1, out);  
     fwrite(&trace.igc, 2, 1, out);  
     fwrite(&trace.igi, 2, 1, out);  
     fwrite(&trace.corr, 2, 1, out); 
     fwrite(&trace.sfs, 2, 1, out);  
     fwrite(&trace.sfe, 2, 1, out);  
     fwrite(&trace.slen, 2, 1, out); 
     fwrite(&trace.styp, 2, 1, out);  
     fwrite(&trace.stas, 2, 1, out);  
     fwrite(&trace.stae, 2, 1, out);  
     fwrite(&trace.tatype, 2, 1, out);  
     fwrite(&trace.afilf, 2, 1, out);  
     fwrite(&trace.afils, 2, 1, out);  
     fwrite(&trace.nofilf, 2, 1, out);  
     fwrite(&trace.nofils, 2, 1, out);  
     fwrite(&trace.lcf, 2, 1, out);  
     fwrite(&trace.hcf, 2, 1, out);  
     fwrite(&trace.lcs, 2, 1, out);  
     fwrite(&trace.hcs, 2, 1, out);  
     fwrite(&trace.year, 2, 1, out);  
     fwrite(&trace.day , 2, 1, out);  
     fwrite(&trace.hour, 2, 1, out);  
     fwrite(&trace.minute, 2, 1, out);
     fwrite(&trace.sec, 2, 1, out);
     fwrite(&trace.timbas, 2, 1, out);
     fwrite(&trace.trwf, 2, 1, out);
     fwrite(&trace.grnors, 2, 1, out);
     fwrite(&trace.grnofr, 2, 1, out);
     fwrite(&trace.grnlof, 2, 1, out);
     fwrite(&trace.gaps, 2, 1, out);
     fwrite(&trace.otrav, 2, 1, out);
     fwrite(&trace.cdpx, 4, 1, out);
     fwrite(&trace.cdpy, 4, 1, out);
     fwrite(&trace.iline, 4, 1, out);
     fwrite(&trace.cline, 4, 1, out);
     fwrite(&trace.shotpt, 4, 1, out);
     fwrite(&messs, 2, 1, out); //201
     fwrite(&messs, 2, 1, out); //203
     fwrite(&mess,  4, 1, out); //205
     fwrite(&messs, 2, 1, out); //209
     fwrite(&messs, 2, 1, out);
     fwrite(&messs, 2, 1, out);
     fwrite(&messs, 2, 1, out);
     fwrite(&messs, 2, 1, out);
     fwrite(&messs, 2, 1, out);
     fwrite(&mess, 4, 1, out);
     fwrite(&mess, 4, 1, out);
     fwrite(&messs, 2, 1, out);
     fwrite(&messs, 2, 1, out);
     fwrite(&mess, 4, 1, out);
     mess=htonl(9999);
     fwrite(&mess, 4, 1, out);

     // writing data to image file
     fwrite(data, 4, ns, out);

   } // if loop over trace number
   it++;
  } // offset loop
   csp=csp+1;
  } //yc
  } // xc
 fclose(out);

//////////////////////////////////////////////////
 free(csp_data); //free_array4d(csp_data,nx,ny,nh,ns);
 free(block_trace);
 free(block_header);

 fclose(input1);
 fclose(inputi);

//////////////////////////////////////////////////
 fftw_destroy_plan( plan_forward );
 fftw_destroy_plan( plan_backward );
 fftw_free( dataT );
 fftw_free( fft_result );
 fftw_free( ifft_result );

}
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////