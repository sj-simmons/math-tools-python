#include <fftw3.h>     /* mpicc thisfile.c -lfftw3 -laudiofile -lm           */
#include <stdio.h>     /* This implements a time-varying filter using filter */
#include <stdlib.h>    /* banks along the lines of, e.g. [1], in parallel.   */ 
#include <math.h>      /* Copyright 2009 Scott Simmons                       */
#include <audiofile.h>
#include <sys/param.h>
#include <string.h>
#include <mpi.h>       /* [1] A Unified Approach to Short-Time Fourier       */
                       /*     Analysis and Synthesis by Jont B. Allen and    */
                       /*     Lawrence R. Rabiner.  Proceeding of IEEE,      */
                       /*     Vol. 65, No.11, 1977.                          */
#define RESULT 1
#define DIETAG 2
#define OMEGA 1024

main(int argc, char **argv)
{
  int rank, numnodes;
  MPI_Init(&argc,&argv);                     /* Initialize mpi */
  MPI_Comm_size(MPI_COMM_WORLD, &numnodes); /* number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);    /* Find out which node running on */ 
  //double chunk[CHUNKSIZE];
  int framechans;
  double *ninbuf, samprate, *outbuf, *cutoffseq;

  if (rank == 0) {   /* On the head node: */
    MPI_Status status;
    char *filename, string[MAXPATHLEN], *nullstr="", *dot, *suffix;
    AFfilehandle infile, outfile;
    AFfilesetup filesetup;
    AFframecount nframes;
    long nchans;
    int i,j,fmt;
    double *outbuf, *retbuf;

    /* Read audio in with libaudiofile */
    if (argc != 2){fprintf(stderr,"Usage:%s <filename>\n",argv[0]);exit(-1);}
    filename = argv[1];
    infile = afOpenFile(filename,"r",AF_NULL_FILESETUP);
    if (!infile){fprintf(stderr,"Failed to open file %s\n",filename);exit(-1);}
    if (dot = strrchr(filename,'.')) {*dot = 0;suffix = dot+1;} else{suffix = nullstr;}
    nframes = afGetFrameCount(infile,AF_DEFAULT_TRACK);
    nchans = afGetChannels(infile,AF_DEFAULT_TRACK);
    fmt = afGetFileFormat(infile,0);
    samprate = afGetRate(infile,AF_DEFAULT_TRACK);
    afSetVirtualSampleFormat(infile, AF_DEFAULT_TRACK, AF_SAMPFMT_DOUBLE, 32);
    printf("Channels: %d, Frames: %d, Format: %d, Sample rate: %f\n",nchans,nframes,fmt,samprate);
    framechans=nframes*nchans;
    //int inbufsize=framechans*sizeof(double); if (!inbuf){fprintf(stderr,"Unable to allocate memory buffer\n");exit(-1);}
    double inbuf[framechans];
    outbuf = malloc(nframes*nchans*sizeof(double)); if (!outbuf){fprintf(stderr,"Unable to allocate memory buffer\n");exit(-1);}
    filesetup = afNewFileSetup();
    afInitFileFormat(filesetup, fmt);
    afInitChannels(filesetup, AF_DEFAULT_TRACK,  nchans); 
    afInitRate(filesetup, AF_DEFAULT_TRACK, samprate);
    sprintf(string, "%s_filtered.%s",filename,suffix);
    printf("creating %s\n",string);
    outfile = afOpenFile(string,"w",filesetup); if (!outfile){fprintf(stderr,"Open failed for %s\n",string);exit(-1);}
    afSetVirtualSampleFormat(outfile, AF_DEFAULT_TRACK, AF_SAMPFMT_DOUBLE, 32);
    while (1) { /* Read samples and count size. */
      int nread; nread = afReadFrames(infile, AF_DEFAULT_TRACK, inbuf, nframes);
      if (!nread) { break; }
    }

    printf("Head broadcasting frames*channels: %d and sample rate: %f:\n",framechans,samprate);
    MPI_Bcast(&framechans,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&samprate,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    printf("Head broadcasting inbuf with size %d:\n",framechans);
    MPI_Bcast(inbuf,framechans,MPI_DOUBLE,0,MPI_COMM_WORLD);

    outbuf = malloc(framechans*sizeof(double));
    retbuf = malloc(framechans*sizeof(double));

    double varcutoff[framechans];

    fftw_complex *signal;
    fftw_plan p;
    int rectwin=2048,k;
    double max,mag;

    signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*rectwin);
    for (k=0;k < floor(framechans/rectwin);k++){
      for (i=0;i<rectwin;i++) { signal[i][1]=inbuf[k*rectwin+i]; signal[i][2]=0.0; }
      p=fftw_plan_dft_1d(rectwin,signal,signal,FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);

      mag = 0.0;
      for (i=0;i<rectwin;i++){
        mag=mag+sqrt(signal[i][1]*signal[i][1]+signal[i][2]*signal[i][2]);
        for (j=0;j<rectwin;j++){
          *(varcutoff+k*rectwin+j)=mag;
        }
      }
    }

    for (i=floor(framechans/rectwin)*rectwin;i<framechans;i++)
      *(varcutoff+i)=500.0;

    for (i=0;i<framechans;i++){
      printf("%d %f..",i,*(varcutoff+i));
    }

    double sum; 
    int blah=1000;
    for (i=0;i<framechans-blah;i++) {
      sum=0.0;
      for (j=0;j<blah;j++) {
        sum = sum + *(varcutoff+i+j);
      }
      *(varcutoff+i)=sum/blah;
    }

    max=0.0;
    double min=20000.0;
    for (i=0;i<framechans;i++){
      if (*(varcutoff+i) > max) { max = *(varcutoff+i);}
      if (*(varcutoff+i) < min) { min = *(varcutoff+i);}
    }

    double temp;
    for (i=0;i<framechans;i++){
       temp = 200.0 + (*(varcutoff+i)-min)*(5000.0-200.0)/(max-min);
       *(varcutoff+i)=temp;
    } 

    /* generate sinusoidal cutoff sequence */
    // for(i=0;i<framechans;i++)
    // *(varcutoff+i)=1200.0-900.0*cos(i/1500.0);  

    printf("Head broadcasting varcutoff with size %d:\n",framechans);
    MPI_Bcast(varcutoff,framechans,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for (i=0;i<framechans;i++){*(outbuf+i)=0.0;}

    for (i=1;i<numnodes;i++)
       MPI_Send(&i,1,MPI_INT,i,RESULT,MPI_COMM_WORLD);

    for (i=numnodes;i<framechans;i++){
       printf("node %d receiving retbuf..",rank);
       MPI_Recv(retbuf,framechans,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
       printf("node %d done receiving retbuf from node %d\n",rank,status.MPI_SOURCE);
       
       for (j=0;j<framechans;j++)
         *(outbuf+i) =+ *(retbuf+i);
       MPI_Send(&i,1,MPI_INT,status.MPI_SOURCE,RESULT,MPI_COMM_WORLD);
    }

    for (i=1;i<numnodes;i++){
      MPI_Recv(retbuf,framechans,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      for (j=0;j<framechans;j++)
        *(outbuf+i) =+ *(retbuf+i);
    }

    for (i=1;i<numnodes;i++)
      MPI_Send(0,0,MPI_INT,i,DIETAG,MPI_COMM_WORLD);

    printf("  writing audiofile out buffer ... ");
    for (i=0;i<framechans;i++) { //amplify
      *(outbuf+i)=*(outbuf+i)*6.0; 
//    *(outbuf+i)=(sqrt(inverseout[i][1]*inverseout[i][1]+inverseout[i][2]*inverseout[i][2])/((double)paddedlen));
    }
    printf("done.\n");

    printf("writing %s\n",string);
    /* write samples to new wave file */
    afWriteFrames(outfile, AF_DEFAULT_TRACK,outbuf, nframes);
    afCloseFile(infile);
    afCloseFile(outfile);

    free(inbuf);
    free(outbuf);
    free(retbuf);
    fftw_free(signal);
  }                /* end head node */

  if (rank > 0) {  /* begin slaves */
    double *hammbuf, *sincbuf, *kernbuf, *retbuf;
    int i,pow2,len,L=1;  /* L = dft size (power of 2) */
    fftw_complex *signal, *kernel, *inverse;
    fftw_plan p;
    MPI_Status status;
    int varcut;
    int winwidth=1024;

    printf("On node %d\n",rank);

    MPI_Bcast(&framechans,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&samprate,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    ninbuf = malloc(framechans*sizeof(double));
    MPI_Bcast(ninbuf,framechans,MPI_DOUBLE,0,MPI_COMM_WORLD);
    cutoffseq = malloc(framechans*sizeof(double));
    MPI_Bcast(cutoffseq,framechans,MPI_DOUBLE,0,MPI_COMM_WORLD);

   // for (i=0;i<100;i++) 
    //  printf("%f ", *(ninbuf+i));

  //  pow2=ceil(log(framechans)/log(2));
  //  for (i=0;i<pow2;i++) {L*=2;}
  //  printf("node %d: next power of 2 = %d\n padding=%d\n",rank,L,L-framechans);
    L = framechans+winwidth;

    if ((framechans%2)==0) {len=framechans+1;} else {len=framechans;}

    void hamm (double *hammp, int windowlen, int len, int L)
    { int i; 
      for (i=0;i<L;i++)
        *(hammp+i)=0.0;
      for (i=0;i<windowlen;i++) 
        *(hammp+(len-windowlen)/2+1+i)=0.53836-0.46164*cos(2.0*M_PI*i/len);
    }

    hammbuf = malloc(L*sizeof(double)); if (!hammbuf) {fprintf(stderr,"Unable to allocate hammbuf memory buffer\n");exit(-1);}
    hamm(hammbuf,winwidth-1,len,L); //use odd second parameter

    sincbuf=malloc(L*sizeof(double)); if (!sincbuf) {fprintf(stderr,"Unable to allocate sincbuf memory buffer\n");exit(-1);}
    retbuf =malloc(L*sizeof(double)); if (!retbuf) {fprintf(stderr,"Unable to allocate retbuf memory buffer\n");exit(-1);}
    kernbuf = malloc(L*sizeof(double)); if (!kernbuf) {fprintf(stderr,"Unable to allocate kernbuf memory buffer\n");exit(-1);}

    signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*L);
    kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*L);
    inverse =(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*L);

    while (1){

      MPI_Recv(&varcut,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      printf("node %d received varcut %d\n",rank,varcut);

      if (status.MPI_TAG == DIETAG) { return;}

      void sinc (double *sincp, int len, int L, double d) /* use odd len */ 
      { int i, shift=(len-1)/2;
        for (i=0;i<len;i++) {
          if (i-shift != 0)
            *(sincp+i)=sin(M_PI*2.0*d*(i-shift))/(M_PI*(i-shift)); 
          else
            *(sincp+i)=2.0*d;
        }
        for (i=len;i<L;i++)
          *(sincp+i)=0.0; 
      }

      sinc(sincbuf,len,L,*(cutoffseq+varcut)/samprate);

      for (i=0;i<L;i++)
        *(kernbuf+i)=(*(hammbuf+i))*(*(sincbuf+i));

      for (i=(len-1)/2;i<L;i++) /* shift kernbuf here */
        kernel[i-(len-1)/2][1]=*(kernbuf+i);
      for (i=0;i<(len-1)/2;i++) 
        kernel[L-(len-1)/2+i][1]=*(kernbuf+i);
      for (i=0;i<L;i++) 
        kernel[i][2]=(double)0.0;

      /* pad signal */
      for (i=0;i<framechans;i++) { signal[i][1]=*(ninbuf+i); signal[i][2]=0.0; }
      for (i=framechans;i<L;i++) { signal[i][1]=0.0; signal[i][2]=0.0; }

      p=fftw_plan_dft_1d(L,signal,signal,FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);

      p=fftw_plan_dft_1d(L,kernel,kernel,FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);

      for (i=0;i<L;i++) {
        inverse[i][1]=(signal[i][1]*kernel[i][1]-signal[i][2]*kernel[i][2]);
        inverse[i][2]=(signal[i][1]*kernel[i][2]+signal[i][2]*kernel[i][1]);
      }

      p=fftw_plan_dft_1d(L,inverse,inverse,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);

      for (i=0;i<framechans;i++)
        *(retbuf+i)=inverse[i][1]/L;

      printf("node %d sending retbuf...",rank);
      MPI_Send(retbuf,framechans,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      printf("node %d done sending retbuf\n",rank);
    }

    printf("  freeing memory ... ");
    fftw_free(inverse); fftw_free(kernel); fftw_free(signal);
    free(sincbuf); free(kernbuf); free(hammbuf);
    printf("done.\n");
  } /* end slaves */

  MPI_Barrier(MPI_COMM_WORLD);
  printf("rank %d finalizing MPI and exiting.\n",rank);
  MPI_Finalize();
}
