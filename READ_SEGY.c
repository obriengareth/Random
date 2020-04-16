//  Gareth O'Brien 12/6/2012
// -- Include intrinsic C calls --
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<arpa/inet.h>
#include"mpi.h"

uint32_t htonl(uint32_t hostlong);
uint16_t htons(uint16_t hostshort);
uint32_t ntohl(uint32_t netlong);
uint16_t ntohs(uint16_t netshort);

float ntohf( const float inFloat );
void count_size_file(char *argv[]);
void read_segy_header(char *argv[], int tracenum);
void read_segy_trace(char *argv[], int tracenum);

/////////////////////////////////////////////////////////////////////////////
int ns,nt,dt;
int my_rank,my_size;
int rank,size;
int dims[1],periods[1],reorder;
const  int  ndims      = 1;
MPI_Comm new_comm;
/////////////////////////////////////////////////////////////////////////////

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


//  MAIN CONTROL LOOP
main(int argc ,char *argv[])
{
 int i,j,k;
 double start,finish;
 char command[50];

// Starting MPI calls and redefining new communicator
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

// command line arguements
  if( argc != 2 )
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
    printf("Reading Segy\n");
    printf("--------------------------------------------------------\n");
  }

  count_size_file(argv);
  read_segy_header(argv,1);
  read_segy_trace(argv,1);   // ieee format only


 finish=MPI_Wtime();
// printf("finished on node %d after %lf minutes\n",my_rank,(finish-start)/60);
 MPI_Barrier(new_comm);
 MPI_Finalize();

}          // End of main


void count_size_file(char *argv[])
{
  unsigned long length;
  FILE *input;
  char file_name[50];
  ushort ns1,dt1;

  sprintf(file_name,"%s",argv[1]);
  if(( input = fopen(file_name,"rb")) == NULL )
  {
    printf("Error 12345: opening %s on %d\n",my_rank);
    MPI_Finalize();
    exit(1);
  }

 fseek(input, 0, SEEK_END);
 length=ftell(input);
 fseek(input, 3200 + 400 + 114, SEEK_SET);

 fread(&ns1 , sizeof(ushort), 1, input);  ns1=htons(ns1);
 fread(&dt1 , sizeof(ushort), 1, input);  dt1=htons(dt1);
 fclose(input);

 ns= (int) ns1;
 dt= (int) dt1;
 nt= ((length-3600)/(240+ns1*4));

 if(my_rank==0)
 {
  //printf("Number of samples = %lu\n" ,ns1);
  //printf("Time step in ms = %lu\n" ,dt1);
  //printf("Size of the file = %lu\n" ,length);
  //printf("Number of traces = %lu\n" , (length-3600)/(240+ns1*4));
  printf("Size of the file = %lu bytes\n" ,length);
  printf("Number of samples = %d\n" ,ns);
  printf("Time step in ms = %d\n" ,dt);
  printf("Number of traces = %d\n" , nt);
  printf("--------------------------------------------------------\n");
 }

}


void read_segy_header(char *argv[], int tracenum)
{

  FILE *input;
  char file_name[50];
  headers trace;

  sprintf(file_name,"%s",argv[1]);
  if(( input = fopen(file_name,"rb")) == NULL )
  {
    printf("Error 12345: opening %s on %d\n",my_rank);
    MPI_Finalize();
    exit(1);
  }

 fseek(input, 3200 + 400  + 240*(tracenum-1)+(tracenum-1)*ns*4, SEEK_SET);

 fread(&trace.tracl, sizeof(int), 1, input);  trace.tracl=htonl(trace.tracl);
 fread(&trace.tracr, sizeof(int), 1, input);  trace.tracr=htonl(trace.tracr);
 fread(&trace.fldr,  sizeof(int), 1, input);  trace.fldr=htonl(trace.fldr);
 fread(&trace.tracf, sizeof(int), 1, input);  trace.tracf=htonl(trace.tracf);
 fread(&trace.ep,    sizeof(int), 1, input);  trace.ep=htonl(trace.ep);
 fread(&trace.cdp,   sizeof(int), 1, input);  trace.cdp=htonl(trace.cdp);
 fread(&trace.cdpt,  sizeof(int), 1, input);  trace.cdpt=htonl(trace.cdpt);

 fread(&trace.trid, sizeof(short), 1, input);  trace.trid=htons(trace.trid);
 fread(&trace.nvs,  sizeof(short), 1, input);  trace.nvs=htons(trace.nvs);
 fread(&trace.nhs,  sizeof(short), 1, input);  trace.nhs=htons(trace.nhs);
 fread(&trace.duse, sizeof(short), 1, input);  trace.duse=htons(trace.duse);

 fread(&trace.offset, sizeof(int), 1, input);  trace.offset=htonl(trace.offset);
 fread(&trace.gelev,  sizeof(int), 1, input);  trace.gelev=htonl(trace.gelev);
 fread(&trace.selev,  sizeof(int), 1, input);  trace.selev=htonl(trace.selev);
 fread(&trace.sdepth, sizeof(int), 1, input);  trace.sdepth=htonl(trace.sdepth);
 fread(&trace.gdel,  sizeof(int), 1, input);   trace.gdel=htonl(trace.gdel);
 fread(&trace.sdel,  sizeof(int), 1, input);   trace.sdel=htonl(trace.sdel);
 fread(&trace.swdep, sizeof(int), 1, input);   trace.swdep=htonl(trace.swdep);
 fread(&trace.gwdep, sizeof(int), 1, input);   trace.gwdep=htonl(trace.gwdep);

 fread(&trace.scalel, sizeof(short), 1, input);  trace.scalel=htons(trace.scalel);
 fread(&trace.scalco, sizeof(short), 1, input);  trace.scalco=htons(trace.scalco);

 fread(&trace.sx, sizeof(int), 1, input);  trace.sx=htonl(trace.sx);
 fread(&trace.sy, sizeof(int), 1, input);  trace.sy=htonl(trace.sy);
 fread(&trace.gx, sizeof(int), 1, input);  trace.gx=htonl(trace.gx);
 fread(&trace.gy, sizeof(int), 1, input);  trace.gy=htonl(trace.gy);

 fread(&trace.counit, sizeof(short), 1, input);  trace.counit=htons(trace.counit);
 fread(&trace.wevel,  sizeof(short), 1, input);  trace.wevel=htons(trace.wevel);
 fread(&trace.swevel, sizeof(short), 1, input);  trace.swevel=htons(trace.swevel);
 fread(&trace.sut,    sizeof(short), 1, input);  trace.sut=htons(trace.sut);
 fread(&trace.gut,    sizeof(short), 1, input);  trace.gut=htons(trace.gut);
 fread(&trace.sstat,  sizeof(short), 1, input);  trace.sstat=htons(trace.sstat);
 fread(&trace.gstat,  sizeof(short), 1, input);  trace.gstat=htons(trace.gstat);
 fread(&trace.tstat,  sizeof(short), 1, input);  trace.tstat=htons(trace.tstat);
 fread(&trace.laga,   sizeof(short), 1, input);  trace.laga=htons(trace.laga);
 fread(&trace.lagb,   sizeof(short), 1, input);  trace.lagb=htons(trace.lagb);
 fread(&trace.delrt,  sizeof(short), 1, input);  trace.delrt=htons(trace.delrt);
 fread(&trace.muts,   sizeof(short), 1, input);  trace.muts=htons(trace.muts);
 fread(&trace.mute ,  sizeof(short), 1, input);  trace.mute=htons(trace.mute );

 fread(&trace.ns , sizeof(ushort), 1, input);  trace.ns=htons(trace.ns);
 fread(&trace.dt , sizeof(ushort), 1, input);  trace.dt=htons(trace.dt);

 fread(&trace.gain, sizeof(short), 1, input);  trace.gain=htons(trace.gain);
 fread(&trace.igc, sizeof(short), 1, input);  trace.igc=htons(trace.igc);
 fread(&trace.igi, sizeof(short), 1, input);  trace.igi=htons(trace.igi);
 fread(&trace.corr, sizeof(short), 1, input);  trace.corr=htons(trace.corr);
 fread(&trace.sfs, sizeof(short), 1, input);  trace.sfs=htons(trace.sfs);
 fread(&trace.sfe, sizeof(short), 1, input);  trace.sfe=htons(trace.sfe);
 fread(&trace.slen, sizeof(short), 1, input);  trace.slen=htons(trace.slen);
 fread(&trace.styp, sizeof(short), 1, input);  trace.styp=htons(trace.styp);
 fread(&trace.stas, sizeof(short), 1, input);  trace.stas=htons(trace.stas);
 fread(&trace.stae, sizeof(short), 1, input);  trace.stae=htons(trace.stae);
 fread(&trace.tatype, sizeof(short), 1, input);  trace.tatype=htons(trace.tatype);
 fread(&trace.afilf, sizeof(short), 1, input);  trace.afilf=htons(trace.afilf);
 fread(&trace.afils, sizeof(short), 1, input);  trace.afils=htons(trace.afils);
 fread(&trace.nofilf, sizeof(short), 1, input);  trace.nofilf=htons(trace.nofilf);
 fread(&trace.nofils, sizeof(short), 1, input);  trace.nofils=htons(trace.nofils);
 fread(&trace.lcf, sizeof(short), 1, input);  trace.lcf=htons(trace.lcf);
 fread(&trace.hcf, sizeof(short), 1, input);  trace.hcf=htons(trace.hcf);
 fread(&trace.lcs, sizeof(short), 1, input);  trace.lcs=htons(trace.lcs);
 fread(&trace.hcs, sizeof(short), 1, input);  trace.hcs=htons(trace.hcs);
 fread(&trace.year, sizeof(short), 1, input);  trace.year=htons(trace.year);
 fread(&trace.day , sizeof(short), 1, input);  trace.day =htons(trace.day );
 fread(&trace.hour, sizeof(short), 1, input);  trace.hour=htons(trace.hour);
 fread(&trace.minute, sizeof(short), 1, input);  trace.minute=htons(trace.minute);
 fread(&trace.sec, sizeof(short), 1, input);  trace.sec=htons(trace.sec);
 fread(&trace.timbas, sizeof(short), 1, input);  trace.timbas=htons(trace.timbas);
 fread(&trace.trwf, sizeof(short), 1, input);  trace.trwf=htons(trace.trwf);
 fread(&trace.grnors, sizeof(short), 1, input);  trace.grnors=htons(trace.grnors);
 fread(&trace.grnofr, sizeof(short), 1, input);  trace.grnofr=htons(trace.grnofr);
 fread(&trace.grnlof, sizeof(short), 1, input);  trace.grnlof=htons(trace.grnlof);
 fread(&trace.gaps, sizeof(short), 1, input);  trace.gaps=htons(trace.gaps);
 fread(&trace.otrav, sizeof(short), 1, input);  trace.otrav=htons(trace.otrav);
 //181 - 240     Unassigned for optional info
 fread(&trace.cdpx, sizeof(int), 1, input);  trace.cdpx=htonl(trace.cdpx);
 fread(&trace.cdpy, sizeof(int), 1, input);  trace.cdpy=htonl(trace.cdpy);
 fread(&trace.iline, sizeof(int), 1, input);  trace.iline=htonl(trace.iline);
 fread(&trace.cline, sizeof(int), 1, input);  trace.cline=htonl(trace.cline);
 fread(&trace.shotpt, sizeof(int), 1, input);  trace.shotpt=htonl(trace.shotpt);
//201 free for moment

 fclose(input);

/* if(my_rank==0)
 {
  printf("%d\n" ,trace.tracl);
  printf("%d\n" ,trace.tracr);
  printf("%d\n" ,trace.fldr);
  printf("%d\n" ,trace.tracf);
  printf("%d\n" ,trace.ep);
  printf("%d\n" ,trace.cdp);
  printf("%d\n" ,trace.cdpt);

  printf("%lu\n" ,trace.trid);
  printf("%lu\n" ,trace.nvs);
  printf("%lu\n" ,trace.nhs);
  printf("%lu\n" ,trace.duse);

  printf("%d\n" ,trace.offset);
  printf("%d\n" ,trace.gelev);
  printf("%d\n" ,trace.selev);
  printf("%d\n" ,trace.sdepth);
  printf("%d\n" ,trace.gdel);
  printf("%d\n" ,trace.sdel);
  printf("%d\n" ,trace.swdep);
  printf("%d\n" ,trace.gwdep);
  
  printf("%lu\n" ,trace.scalel);
  printf("%lu\n" ,trace.scalco);

  printf("%d\n" ,trace.sx);
  printf("%d\n" ,trace.sy);
  printf("%d\n" ,trace.gx);
  printf("%d\n" ,trace.gy);

  printf("%lu\n" ,trace.counit);
  printf("%lu\n" ,trace.wevel);
  printf("%lu\n" ,trace.swevel);
  printf("%lu\n" ,trace.sut);
  printf("%lu\n" ,trace.gut);
  printf("%lu\n" ,trace.sstat);
  printf("%lu\n" ,trace.gstat);
  printf("%lu\n" ,trace.tstat);
  printf("%lu\n" ,trace.laga);
  printf("%lu\n" ,trace.lagb);
  printf("%lu\n" ,trace.delrt);
  printf("%lu\n" ,trace.muts);
  printf("%lu\n" ,trace.mute);

  printf("%lu\n" ,trace.ns);
  printf("%lu\n" ,trace.dt);

  printf("%lu\n" ,trace.gain);
  printf("%lu\n" ,trace.igc);
  printf("%lu\n" ,trace.igi);
  printf("%lu\n" ,trace.corr);
  printf("%lu\n" ,trace.sfs);
  printf("%lu\n" ,trace.sfe);
  printf("%lu\n" ,trace.slen);
  printf("%lu\n" ,trace.styp);
  printf("%lu\n" ,trace.stas);
  printf("%lu\n" ,trace.stae);
  printf("%lu\n" ,trace.tatype);
  printf("%lu\n" ,trace.afilf);
  printf("%lu\n" ,trace.afils);
  printf("%lu\n" ,trace.nofilf);
  printf("%lu\n" ,trace.nofils);
  printf("%lu\n" ,trace.lcf);
  printf("%lu\n" ,trace.hcf);
  printf("%lu\n" ,trace.lcs);
  printf("%lu\n" ,trace.hcs);
  printf("%lu\n" ,trace.year);
  printf("%lu\n" ,trace.day);
  printf("%lu\n" ,trace.hour);
  printf("%lu\n" ,trace.minute);
  printf("%lu\n" ,trace.sec);
  printf("%lu\n" ,trace.timbas);
  printf("%lu\n" ,trace.trwf);
  printf("%lu\n" ,trace.grnors);
  printf("%lu\n" ,trace.grnofr);
  printf("%lu\n" ,trace.grnlof);
  printf("%lu\n" ,trace.gaps);
  printf("%lu\n" ,trace.otrav);

  printf("%d\n" ,trace.cdpx);
  printf("%d\n" ,trace.cdpy);
  printf("%d\n" ,trace.iline);
  printf("%d\n" ,trace.cline);
  printf("%d\n" ,trace.shotpt);
 }
*/


}


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

void read_segy_trace(char *argv[], int tracenum)
{
  int i;
  float data[ns];
  FILE *input;
  char file_name[50];
  headers trace;

  sprintf(file_name,"%s",argv[1]);
  if(( input = fopen(file_name,"rb")) == NULL )
  {
    printf("Error 12345: opening %s on %d\n",my_rank);
    MPI_Finalize();
    exit(1);
  }

 fseek(input, 3600 + 240 + 240*(tracenum-1)+(tracenum-1)*ns*4, SEEK_SET);
 fread(data, sizeof(float), ns, input);

 /*if(my_rank==0)
 {
  for(i=0; i<ns; i++)
  printf(" %f\n" ,ntohf(data[i]));
 }*/
 fclose(input);
}


/* 
/// write individual trace header
     i = 1001;     s = 1;
     fwrite(&n, 4, 1, outX); //0
     fwrite(&n, 4, 1, outX); //4
     fwrite(&i, 4, 1, outX); //8
     fwrite(&n, 4, 1, outX); //12
     fwrite(&dummy1, 12, 1, outX); //16
     fwrite(&s, 2, 1, outX); //28
     fwrite(&dummy1, 6, 1, outX); //30
     fwrite(&n, 4, 1, outX); //36
     fwrite(&dummy1, 74, 1, outX); //40
     s = (short)Max_Time;
     fwrite(&s, 2, 1, outX); //114
     s = (short)(dt*1000000.0);
     fwrite(&s, 2, 1, outX); //116
     fwrite(&dummy1, 122, 1, outX);//118 - 240
*/
