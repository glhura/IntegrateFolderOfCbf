#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>

#include <sys/time.h>
#include <sys/resource.h>

double get_time();
double starttime;
int append_output = 0;

/* POINTER defines a generic pointer type */
typedef unsigned char *POINTER;

/* UINT2 defines a two byte word */
typedef unsigned short int UINT2;

/* UINT4 defines a four byte word */
typedef unsigned int UINT4;

/* MD5 context. */
typedef struct {
  UINT4 state[4];                                   /* state (ABCD) */
  UINT4 count[2];        /* number of bits, modulo 2^64 (lsb first) */
  unsigned char buffer[64];                         /* input buffer */
} MD5_CTX;

void MD5Init(MD5_CTX *);
void MD5Update(MD5_CTX *, unsigned char *, unsigned int);
void MD5Final(unsigned char [16], MD5_CTX *);
int cbf_md5digest_to64 (char *, const unsigned char *);
double CriterionPhi(double x, double A, double B, double C, double a, double b , double c);
double Polarization(double theta , double phi);
double ThicknessCorrection(double theta , double phi, double thickness, double AbsoCoef, double psi);
double AirAbsorption(double r, double PixelSize, double AirAbso);
double PixelThickness (double r , double theta , double phi, double detx, double dety, double detz, double c);
double DotProduct (double r , double theta , double phi, double detx, double dety, double detz);



double **intensity_matrix_gen(char *InFileName, size_t num_of_params, double *expt_params, size_t lines_in_mask, int** mask)
{
    int testinc1 = 0;
    int testinc2 = 0;
    double thetamax = 0.1;
    double thetaResol  = 0.0001;
    double PixelSize = 0.0172;
    int skipone = 0;
    int n,i,j,k,pixels;
    int header_lines=0, outheader_lines=0;
    int header_bytes=0, outheader_bytes=512;
    int file_bytes_read=0, data_bytes_read=0;
    int compressed_size=0, uncompressed_size=0;
    int8_t *inimage1;
    int32_t *inimage2;
    int32_t *decompressed_data;
    char *start_of_data;
    int           PrintPhi_flag = 0;
    int print_timing = 0;
    int overflows=0;
    int underflows=0;
    int bits = 32;
    int bytes = 4;
    int sign = 1;
    int floatout = 0;
    int downsample = 0;
    int nonlinear = 0;
    int nomd5check = 0;
    int md5checkonly = 0;
    int nowrite = 0;
    int nostats = 0;
    float scale=1, offset=0;
    FILE          *infile1;
    char          FileRoot[100] = "";
    char          outfilename[100] = "";
    FILE          *outfile1;
    union {int8_t *int8ptr; uint8_t *uint8ptr;
          int16_t *int16ptr; uint16_t *uint16ptr;
	        int32_t *int32ptr; float *floatptr;
          double *doubleptr;} bitshuffle;
    int32_t cval,cdiff;
    char *headerstuff, *string;
    int max_x=0,max_y=0,xpixels,ypixels;
    float sum,sumd,sumsq,sumdsq,avg,diff,rms,rmsd;
    long max,min,isum=0,valid_pixels=0;
    double        **IntenFromDet;
    int xCompany = 0, yCompany = 0;
    int           StartCorner,ApplyThickness;
    double        xc, yc, dist, TiltPlane, RotAng;
    double        lambda, qmin_set, qmax_set, AbsoCoef, thickness;
    double	      detinfo1, detinfo2, detinfo3, detinfo4;
    double        A, B, C, a, b, c;
    double        ThetaCirc, Q;
    double        r, phi, toppart, ki;
    double        cb, sb, ca;
    double        temp,top,bot,mid,detx,dety,detz,norm,DetR;
    double        Degrees[10], Radians[10], AcceptRad = 0.0005;
    unsigned int           SafetyInc = (round(thetamax/thetaResol));
    double        *StandDevTheta, *AveIntenThetaBin, *CleanInten, *IntenToReport, *ErrorToReport;
    int		        *NumPixels1, *NumPixels2;
    double        **theta;
    FILE          *phi0173, *phi02, *phi03, *phi04, *phi05, *phi08, *phi11, *phi14, *phi17, *phi19;
    int		        inc_mask= 0;
    int		        mask_x = mask[inc_mask][0];
    int		        mask_y = mask[inc_mask][1];
    double	      mask_pixels = 100000000;
    int           x, y, xpr, ypr, zpr, xprpr, yprpr, zprpr;
    int		        thetabin = 0;
    double        Check;

    infile1 = fopen(InFileName, "r" );
   
    if ( !infile1 ) 
    {
        perror( "Error opening ALS data file" );
        exit(1);
    } 
    /*printf("%s\n",InFileName);*/
  
    strncpy(FileRoot,InFileName,((strlen(InFileName)-4)));
    sprintf(outfilename,"%s%s", FileRoot, ".dat");
    outfile1 = fopen (outfilename, "w");

  StartCorner = (int) expt_params[0]; ApplyThickness = (int) expt_params[1];
  xc = expt_params[2]; yc = expt_params[3]; dist = expt_params[4]; TiltPlane = expt_params[5]; RotAng = expt_params[6];
  lambda = expt_params[7]; qmin_set = expt_params[8]; qmax_set = expt_params[9]; AbsoCoef = expt_params[10]; thickness = expt_params[11];
  detinfo1 = expt_params[12]; detinfo2 = expt_params[13]; detinfo3 = expt_params[14]; detinfo4 = expt_params[15];
  TiltPlane = TiltPlane*3.141592/180;

starttime = get_time();

/* follow-up calcs */
if(floatout && bits<32) bits = 32;
bytes = bits/8;
if(downsample == 0) downsample = 1;

/* start reading the cbf image */

/* load what should be more than a header on the first try */
headerstuff = calloc(sizeof(char),65535);
file_bytes_read = fread(headerstuff,sizeof(char),65534,infile1);
//printf("GOTHERE: about to overwrite char:%c%c\n",headerstuff[file_bytes_read-1],headerstuff[file_bytes_read]);
/* make sure this is a zero-terminated string? */
headerstuff[file_bytes_read] = 0;
//printf("GOTHERE: read in %d %d-byte blocks\n",file_bytes_read,sizeof(char));
string=strstr(headerstuff,"\f\032\004\325");
while( string == NULL )
    {
        /* TBD: keep reading if we don't see it? */
        printf("This does not look like a CBF file %s\n",string);
        exit(9);
    }
  /* how many bytes are in the header? */
  header_bytes = string - headerstuff + 4;
  start_of_data = string + 4;
  /* calculate number of data bytes already read in as part of the "header" */
  data_bytes_read = file_bytes_read - header_bytes;
  //printf("GOTHERE: start of data chars:%c%c%c%c\n",start_of_data[0],start_of_data[1],start_of_data[2],start_of_data[3]);

  /* now use the header to figure out how many bytes we still need to read */
  string = strstr(headerstuff,"X-Binary-Size:");
  if(string == NULL)
    {
	  printf("*** Bad header - cannot find X-Binary-Size\n");
    exit(9);
    }
  string += strlen("X-Binary-Size:");
  sscanf(string, "%d", &compressed_size);
  /*printf(" compressed_size = %d\n",compressed_size);*/
  inimage1 = calloc(compressed_size,sizeof(char));

  /* copy over stuff we already read up to the current point in the input file */
  memcpy(inimage1,start_of_data,data_bytes_read);
  /* now finish the read */
  fread(inimage1+data_bytes_read,sizeof(char),compressed_size-data_bytes_read,infile1);
  /* and that is all we need to know */
  fclose(infile1);

  /* how big is the output file going to be? */
  string = strstr(headerstuff,"X-Binary-Number-of-Elements:");
  if(string != NULL)
    {
	  string += strlen("X-Binary-Number-of-Elements:");
	  sscanf(string, "%d", &pixels);
    uncompressed_size = pixels*sizeof(int32_t);
    /*printf(" uncompressed_size = %d\n",uncompressed_size);*/
    }

  //    if(pixels == 0) pixels = 6224001;
  inimage2 = calloc(pixels,sizeof(int32_t));

  /* x-y size */
  string = strstr(headerstuff,"X-Binary-Size-Fastest-Dimension:");
  if(string != NULL)
    {
	  string += strlen("X-Binary-Size-Fastest-Dimension:");
	  sscanf(string, "%d", &xpixels);
    /*printf(" x size = %d\n",xpixels);*/
    }
  string = strstr(headerstuff,"X-Binary-Size-Second-Dimension:");
  if(string != NULL)
    {
	  string += strlen("X-Binary-Size-Second-Dimension:");
	  sscanf(string, "%d", &ypixels);
    /*printf(" y size = %d\n",ypixels);*/
    }
  
  IntenFromDet = (double **) calloc(xpixels, sizeof(double*));
  for(j=0; j< xpixels; j++){
    IntenFromDet[j] =  (double *) calloc(ypixels, sizeof(double));}


	MD5_CTX context;
	unsigned char rawdigest[17];
	char digest[25];
	char header_digest[25];
  string += strlen("Content-MD5: ");

  sum = isum = 0;
  max=0;
  i = 0;
  min=INT_MAX;
  /*printf("decompressing...\n");*/
  /* intialize pointer union to start of compressed data */
  bitshuffle.int8ptr = inimage1;
  decompressed_data = inimage2;
  cval = 0;
 while (bitshuffle.int8ptr-inimage1 < compressed_size)
    {
	  /*
	  ** If the file is corrupted with a lot of zeroes, 
	  ** then (*t.ucp != 0x80) is true too many times,
	  ** and 'dp' tries to go beyond the buffer end.
    so: if ((char *)dp >= (ObjPtr->Data+size))
	  */
	  if (decompressed_data >= inimage2+uncompressed_size) break;
	
	  if (*bitshuffle.uint8ptr != 0x80)
      {
      cdiff = (int)*bitshuffle.int8ptr++;
      }
	  else
	  {

		  bitshuffle.int8ptr++;
		  if (*bitshuffle.int16ptr==0)
		    {
  
		    cval = 0;
  
		    bitshuffle.int16ptr++;
  
		    continue;
		    }
		  if (*bitshuffle.uint16ptr != 0x8000)
        {
  
        cdiff = (int)*bitshuffle.int16ptr++;
        }
		  else
        {
 
        bitshuffle.int16ptr++;

        cdiff = *bitshuffle.int32ptr++;
		    }
      }
  /*printf("Here are some values %d %d\n", cval, cdiff);*/
	/*cval = cval + cdiff;*/
  cval = *decompressed_data++ = cval + cdiff;
  /*if(i < 10){printf("%d\n", cval);}*/
  /*printf("Here we are %d %d %d\n", xCompany, yCompany, cval);*/

  switch(StartCorner)
	  {
	   					case 1: {x = xCompany; y = yCompany; break;}
	   					case 2: {x = (xpixels - 1)  - xCompany; y = yCompany; break;}
	   					case 3: {x = xCompany; y = (ypixels - 1) - yCompany; break;}
	   					case 4: {x = (xpixels - 1) - xCompany; y = (ypixels - 1) - yCompany; break;}
	   					case 5: {x = yCompany; y = xCompany; break;}
	   					case 6: {x = (ypixels - 1)  - yCompany; y = xCompany; break;}
	   					case 7: {x = yCompany; y = (xpixels - 1) - xCompany; break;}
	   					case 8: {x = (ypixels - 1) - yCompany; y = (xpixels - 1) - xCompany; break;}
	  }

  IntenFromDet[x][y] = (double) cval;
  
  if (xCompany == mask_x && yCompany ==  mask_y)
		{

			IntenFromDet[x][y] = mask_pixels;
				if (inc_mask < lines_in_mask - 1){
							inc_mask = inc_mask + 1;}
						mask_x = mask[inc_mask][0];
						mask_y = mask[inc_mask][1];
					}

   if(xCompany < xpixels - 1)
  {
    xCompany = xCompany + 1;    
  }
  else
  {
    xCompany = 0;  
    yCompany = yCompany + 1;  
  }

  ++i;
  
	}

/*END OF DECOMPRESSION OF CBF FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Strating Integration - Probably should split this function but on 03/06/2024 didnt have time*/


  
  RotAng = RotAng*3.141592/180;

  dist = dist/PixelSize/10;

  a = xc - dist*sin(-TiltPlane)*sin(RotAng);
  b = yc + dist*cos(-TiltPlane)*sin(RotAng);
  c = dist*cos(RotAng);
  A = sin(-TiltPlane)*sin(RotAng);
  B = -cos(-TiltPlane)*sin(RotAng);
  C = -cos(RotAng);
  
  cb = -C / sqrt(C*C + B*B);
  sb = B / sqrt(C*C + B*B);
  ca = -(sqrt(1 - A*A));
  temp = sqrt(C*C + B*B);
  top = -A*C/temp;
  mid = -B/temp;
  bot = -C;
  detx = top;
  dety = mid;
  detz = bot;
  norm = sqrt(detx*detx + dety*dety + detz*detz);
  detx = detx/norm;
  dety = dety/norm;
  detz = detz/norm;
  StandDevTheta = (double*)calloc(SafetyInc,sizeof(double));
  AveIntenThetaBin = (double*)calloc(SafetyInc,sizeof(double));
  CleanInten = (double*)calloc(SafetyInc,sizeof(double));
  IntenToReport = (double*)calloc(SafetyInc,sizeof(double));
  ErrorToReport = (double*)calloc(SafetyInc,sizeof(double));
  NumPixels1 = (int*) calloc(SafetyInc,sizeof(int));
  NumPixels2 = (int*) calloc(SafetyInc,sizeof(int));

  
  for(i = 0; i < SafetyInc; i++)
    {
	StandDevTheta[i] = 0.0;
	AveIntenThetaBin[i] = 0.0;
	CleanInten[i] = 0.0;
	IntenToReport[i] = 0.0;
	ErrorToReport[i] = 0.0;
	NumPixels1[i] = 0;
	NumPixels2[i] = 0;
     }


  theta = (double **) calloc(xpixels, sizeof(double*));
  for(j=0; j< xpixels; j++){
    theta[j] = (double*) calloc(ypixels, sizeof(double));}

  if(PrintPhi_flag == 1)
  {
    phi0173 = fopen("phi0173" , "w"); Degrees[0] = 0.173; 
    phi02 = fopen("phi02" , "w"); Degrees[1] = 0.2;
    phi03 = fopen("phi03" , "w"); Degrees[2] = 0.3;
    phi04 = fopen("phi04" , "w"); Degrees[3] = 0.4;
    phi05 = fopen("phi05" , "w"); Degrees[4] = 0.5;
    phi08 = fopen("phi08" , "w"); Degrees[5] = 0.8;
    phi11 = fopen("phi11" , "w"); Degrees[6] = 1.1;
    phi14 = fopen("phi14" , "w"); Degrees[7] = 1.4;
    phi17 = fopen("phi17" , "w"); Degrees[8] = 1.7;
    phi19 = fopen("phi19" , "w"); Degrees[9] = 1.9;
  }
 

for(y = 0; y < ypixels; y++)
    {
       for(x = 0; x < xpixels; x++)
	      {
	      /*Define r Easily*/
	      r = sqrt((double)(c*c + (a - ((double) x))*(a - ((double) x)) + 
		    (b - ((double) y))*(b - ((double) y))));

/* Coordinate Transformation a translation and rotation from an */
/* origin defined by bottom corner of the detector to an origin defined by */
/* the sample holder and the incident x-ray beam*/

	    xpr = (x - a);
	    ypr = (y - b)*cb - c*sb;
	    zpr = -(y - b)*sb - c*cb;
	    xprpr = xpr*ca - zpr*A;
	    yprpr = ypr;
	    zprpr = xpr*A + zpr*ca;
/* Defining theta*/

	    toppart = (A * ((double) x - a) + B * ((double) y - b) - C * c);
	    ki = toppart/r;
	    theta[x][y] = acos((double) ki);



       
    
/* Define phi but must take into account that acos function only ranges from */
/* phi > 0 to phi < Pi so must create an x-y line on the detectoAgbehanate1.ASCr after which */
/* another rule comes into effect */


	    if(theta[x][y] > 0.000)
	    {
	      if(y > CriterionPhi((double) x, A, B, C, a, b , c))
		      phi = acos(xprpr / r / sin(theta[x][y]));
	      else		
        /*2pi is 6.283185*/
		      phi = 6.283185 -  acos(xprpr / r / sin(theta[x][y]));	      
	    }
	    else phi = 0.0000;



/* This if statement has a large mask in it for things like the beam stop or */
/* bad pixels*/
	    DetR = (double) ypixels/2;

      if (theta[x][y] > 0.000001)

	    {
        /*if(x==728 && y == 805){printf("%lf\n", IntenFromDet[x][y]);}*/
        if(IntenFromDet[x][y] < mask_pixels)
	      {
	        IntenFromDet[x][y] = IntenFromDet[x][y]*pow(r/c,2.0);
        /*if(x == 0){printf("%d %d\n", x, y );}*/
	      }

	    thetabin = (round(theta[x][y]/thetaResol));
      /*if (thetabin == 41)
      {printf("%d %d %lf\n",x,y,IntenFromDet[x][y]);
      testinc2 = testinc2+1;}*/


	if(IntenFromDet[x][y] < mask_pixels  && thetabin < SafetyInc)
		{
		  AveIntenThetaBin[thetabin] += IntenFromDet[x][y];
		  NumPixels1[thetabin]++;
      testinc1 = testinc1 + 1;
		 }
  else
	  {
			AveIntenThetaBin[thetabin] += 0;
			NumPixels1[thetabin] += 0;

		}

	   }


	  if(PrintPhi_flag == 1 && IntenFromDet[x][y] < mask_pixels)
	    {
	      if(theta[x][y] > Radians[0] && theta[x][y] < Radians[0]+AcceptRad)
		fprintf(phi0173, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	      if(theta[x][y] > Radians[1] && theta[x][y] < Radians[1]+AcceptRad)
		fprintf(phi02, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	      if(theta[x][y] > Radians[2] && theta[x][y] < Radians[2]+AcceptRad)
		fprintf(phi03, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	      if(theta[x][y] > Radians[3] && theta[x][y] < Radians[3]+AcceptRad)
		fprintf(phi04, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	      if(theta[x][y] > Radians[4] && theta[x][y] < Radians[4]+AcceptRad)
		fprintf(phi05, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	      if(theta[x][y] > Radians[5] && theta[x][y] < Radians[5]+AcceptRad)
		fprintf(phi08, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	      if(theta[x][y] > Radians[6] && theta[x][y] < Radians[6]+AcceptRad)
		fprintf(phi11, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	      if(theta[x][y] > Radians[7] && theta[x][y] < Radians[7]+AcceptRad)
		fprintf(phi14, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	      if(theta[x][y] > Radians[8] && theta[x][y] < Radians[8]+AcceptRad)
		fprintf(phi17, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	      if(theta[x][y] > Radians[9] && theta[x][y] < Radians[9]+AcceptRad)
		fprintf(phi19, "%lf %lf\n" , phi, IntenFromDet[x][y]); 
	    }
	}
  }
  /*printf("This Just In %d %lf \n", testinc2, AveIntenThetaBin[41] );*/
  


  for(i = 0; i < SafetyInc; i++)
    {
    if(NumPixels1[i] > 0) 
    {
      AveIntenThetaBin[i] = (AveIntenThetaBin[i]/((double) NumPixels1[i]));

     }   
    else 
    	AveIntenThetaBin[i] = 0.0;
    /*printf("This is the Ave %lf\n", AveIntenThetaBin[i]);*/
    }	     
				

 for(y = 0; y < ypixels; y++)
    {

       for(x = 0; x < xpixels; x++)
	{
	  
	  thetabin = (round(theta[x][y]/thetaResol));
	  
	  i= (int)(thetabin);
	  
	    if(IntenFromDet[x][y] < mask_pixels && i < SafetyInc)
	    {  
	        
	      StandDevTheta[i] += pow((IntenFromDet[x][y] - AveIntenThetaBin[i]),2);
	      
	      NumPixels2[i] = NumPixels2[i] + 1;
	    }
	  else
	    {
	      StandDevTheta[i] += 0.0;
	      NumPixels2[i] = NumPixels2[i] + 0;
	     }
	   
	}
    }

 
  
  for(i = 0; i < SafetyInc; i++)
    {
      if(NumPixels2[i] > 5)	
      	{
           StandDevTheta[i] = (pow((StandDevTheta[i]/((double)NumPixels2[i])),0.5));
	   ErrorToReport[i] = StandDevTheta[i]/(pow(NumPixels2[i],0.5));
        }
      else
      	{
          StandDevTheta[i] = 10000000001;
	  ErrorToReport[i] = 10000000001;
	 }

    }					
  
  


  for(i = 0; i < SafetyInc; i++)
    {
       NumPixels1[i] = 0;
    }
    
  
 
  for(y = 0; y < ypixels; y++)
    {
       for(x = 0; x < xpixels; x++)
	{
	  thetabin = (round(theta[x][y]/thetaResol));
	  i= (int)thetabin;
 
	  if(IntenFromDet[x][y] < mask_pixels && NumPixels2[i] > 5 && i < SafetyInc)
	    {
	      
	    Check = pow((AveIntenThetaBin[i] - IntenFromDet[x][y])*(AveIntenThetaBin[i] - IntenFromDet[x][y]),0.5);
	    if(Check < 5*StandDevTheta[i])
	      	{
		        
	     		CleanInten[i] += IntenFromDet[x][y];
		
	     		NumPixels1[i] = (int)(NumPixels1[i] + 1);
	      	}
		
	      }
	    	    }
	}



 for(i = 0; i < SafetyInc; i++)
    {
	
	if(NumPixels1[i] > 0)
	{	
    		IntenToReport[i] = CleanInten[i]/((double) NumPixels1[i]);
	}
	if(NumPixels1[i] <= 0)
	{
		IntenToReport[i] = 0.0;
	}

    }


 for(i = 0; i < SafetyInc; i++)
 	{
  		if(ErrorToReport[i] == 10000000001)
			{
				for(j = 0; j < SafetyInc; j++)
					{
						if(ErrorToReport[j] !=  10000000001)
							{
							    
							     ErrorToReport[i]  = ErrorToReport[j];
							     j = SafetyInc;
							 }}}} 
   
  for(i = 0; i< SafetyInc; i++)
    {
         ThetaCirc = (double) (((double)i+0.5)*(double)thetaResol);
        
	       Q = ((12.56637/lambda)*sin(ThetaCirc/2.0));

	 if(Q > qmin_set && Q < qmax_set)
		{
       fprintf(outfile1, "%lf %lf %lf\n", Q, IntenToReport[i], ErrorToReport[i]);
			}
    }
    
         /*4pi is 12.56637*/
        /*printf("%lf %lf %lf\n", Q, IntenToReport[i], ErrorToReport[i]);*/
        /*if (Q > 0.358 && Q < 0.37)*/
	 	    /*printf("%lf %lf %lf\n", Q, qmin_set, qmax_set);*/

fclose(outfile1);

if(PrintPhi_flag == 1)
{
	fclose(phi0173);
	fclose(phi02);
	fclose(phi03);
	fclose(phi04);
	fclose(phi05);
	fclose(phi08);
	fclose(phi11);
	fclose(phi14);
	fclose(phi17);
	fclose(phi19);
}

return (IntenFromDet);
  
}

int main (int argc, char *argv[])
{
  unsigned int		mask_lines = 0;
  int           StartCorner,ApplyThickness;
  int 		ch=0;
  int		**MaskPixels;
  int 		size_params = 16;
  double    expt_params[size_params];
  int		count;
  double	param_value;
  int		mask_x, mask_y;
  int		j;
  double	qmin_set, qmax_set;
  double    xc, yc, dist, TiltPlane, RotAng;
  double    lambda, AbsoCoef, thickness, psi, AirAbso;
  char      *InFileName;
  FILE		*MaskFile1;
  FILE      *Params;

  if ( argc != 2 ) {
    fprintf( stderr, "Usage: %s Dectrisfile.tif \n" );
    exit(1);
  }
  
  InFileName = argv[1];
    
  Params = fopen("ExptParams", "r");
  if ( !Params ) 
  {
    perror( "Error opening ExptParams data file" );
    exit(1);
  }
  
  fscanf(Params, "%lf %lf", &expt_params[0], &expt_params[1]);
  fscanf(Params, "%lf %lf %lf %lf %lf" , &expt_params[2], &expt_params[3], &expt_params[4], &expt_params[5], &expt_params[6]);
  fscanf(Params, "%lf %lf %lf %lf %lf" , &expt_params[7], &expt_params[8], &expt_params[9], &expt_params[10], &expt_params[11]);
  fscanf(Params, "%lf %lf %lf %lf" , &expt_params[12], &expt_params[13], &expt_params[14], &expt_params[15]);

  MaskFile1 = fopen("X_Y_mask.asc", "r");
  
  if ( !MaskFile1 ) 
  {
    perror( "Error opening X_Y_mask.asc data file" );
    exit(1);
  }

  while ( !feof(MaskFile1))
    {
      ch = fgetc(MaskFile1);
      if (ch == '\n')
      {
    	mask_lines++;
       }
    }

  fclose(MaskFile1);

  MaskPixels = (int**) calloc(mask_lines, sizeof(int*));
  for(j=0; j< mask_lines; j++){
    MaskPixels[j] = (int*) calloc(2, sizeof(int));}
  
  MaskFile1 = fopen("X_Y_mask.asc", "r");
  
  for(j = 0; j < mask_lines; j++)
  	{
		/*printf("%d\n", mask_lines);*/
		fscanf(MaskFile1, "%d %d", &MaskPixels[j][0], &MaskPixels[j][1]);
		/*printf("%d %d\n", MaskData[j][0], MaskData[j][1]);*/

	}
    
 /*Need To figure out how to load data into mask*/
  
  intensity_matrix_gen(InFileName, (size_t) 16, expt_params, (size_t) mask_lines, MaskPixels);
    
 }





/*
Eric F. Eikenberry's MD5 thing that doesn't obey the MD5 standard and nobody semes to use


   MD5.H - header file for MD5C.C
 */

/* Copyright (C) 1991-2, RSA Data Security, Inc. Created 1991. All
rights reserved.

License to copy and use this software is granted provided that it
is identified as the "RSA Data Security, Inc. MD5 Message-Digest
Algorithm" in all material mentioning or referencing this software
or this function.

License is also granted to make and use derivative works provided
that such works are identified as "derived from the RSA Data
Security, Inc. MD5 Message-Digest Algorithm" in all material
mentioning or referencing the derived work.

RSA Data Security, Inc. makes no representations concerning either
the merchantability of this software or the suitability of this
software for any particular purpose. It is provided "as is"
without express or implied warranty of any kind.

These notices must be retained in any copies of any part of this
documentation and/or software.
 */
double get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

/* Constants for MD5Transform routine.
 */
#define S11 7
#define S12 12
#define S13 17
#define S14 22
#define S21 5
#define S22 9
#define S23 14
#define S24 20
#define S31 4
#define S32 11
#define S33 16
#define S34 23
#define S41 6
#define S42 10
#define S43 15
#define S44 21

static void MD5Transform (UINT4 [4], unsigned char [64]);
static void Encode(unsigned char *, UINT4 *, unsigned int);
static void Decode(UINT4 *, unsigned char *, unsigned int);
static void MD5_memcpy(POINTER, POINTER, unsigned int);
static void MD5_memset(POINTER, int, unsigned int);

static unsigned char PADDING[64] = {
  0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/* F, G, H and I are basic MD5 functions.
 */
#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define I(x, y, z) ((y) ^ ((x) | (~z)))

/* ROTATE_LEFT rotates x left n bits.
 */
/* #define ROTATE_LEFT(x, n) (((x) << (n)) | ((x) >> (32-(n)))) */

#define ROTATE_LEFT(x, n) (((x) << (n)) | (((x) & 0x0FFFFFFFF) >> (32 - (n))))

/* FF, GG, HH, and II transformations for rounds 1, 2, 3, and 4.
Rotation is separate from addition to prevent recomputation.
 */
#define FF(a, b, c, d, x, s, ac) { \
 (a) += F ((b), (c), (d)) + (x) + (UINT4)(ac); \
 (a) = ROTATE_LEFT ((a), (s)); \
 (a) += (b); \
  }
#define GG(a, b, c, d, x, s, ac) { \
 (a) += G ((b), (c), (d)) + (x) + (UINT4)(ac); \
 (a) = ROTATE_LEFT ((a), (s)); \
 (a) += (b); \
  }
#define HH(a, b, c, d, x, s, ac) { \
 (a) += H ((b), (c), (d)) + (x) + (UINT4)(ac); \
 (a) = ROTATE_LEFT ((a), (s)); \
 (a) += (b); \
  }
#define II(a, b, c, d, x, s, ac) { \
 (a) += I ((b), (c), (d)) + (x) + (UINT4)(ac); \
 (a) = ROTATE_LEFT ((a), (s)); \
 (a) += (b); \
  }

/* MD5 initialization. Begins an MD5 operation, writing a new context.
 */
void MD5Init (context)
MD5_CTX *context;                                        /* context */
{
  context->count[0] = context->count[1] = 0;
  /* Load magic initialization constants.
*/
  context->state[0] = 0x67452301;
  context->state[1] = 0xefcdab89;
  context->state[2] = 0x98badcfe;
  context->state[3] = 0x10325476;
}

/* MD5 block update operation. Continues an MD5 message-digest
  operation, processing another message block, and updating the
  context.
 */
void MD5Update (context, input, inputLen)
MD5_CTX *context;                                        /* context */
unsigned char *input;                                /* input block */
unsigned int inputLen;                     /* length of input block */
{
  unsigned int i, index, partLen;
  UINT4 I1, I2, S;

  /* Compute number of bytes mod 64 */
  index = (unsigned int)((context->count[0] >> 3) & 0x3F);

  /* Update number of bits */
  I1 = ((UINT4) inputLen) << 3;
  I2 = ((UINT4) context->count [0]);
  context->count[0] = S = I1 + I2;
  if (((~S & (I1 | I2)) | (I1 & I2)) & 0x080000000)
    context->count[1]++;
  context->count[1] += ((UINT4) inputLen >> 29);

  partLen = 64 - index;

  /* Transform as many times as possible.
*/
  if (inputLen >= partLen) {
 MD5_memcpy
   ((POINTER)&context->buffer[index], (POINTER)input, partLen);
 MD5Transform (context->state, context->buffer);

 for (i = partLen; i + 63 < inputLen; i += 64)
   MD5Transform (context->state, &input[i]);

 index = 0;
  }
  else
 i = 0;

  /* Buffer remaining input */
  MD5_memcpy
 ((POINTER)&context->buffer[index], (POINTER)&input[i],
  inputLen-i);
}

/* MD5 finalization. Ends an MD5 message-digest operation, writing the
  the message digest and zeroizing the context.
 */
void MD5Final (digest, context)
unsigned char digest[16];                         /* message digest */
MD5_CTX *context;                                       /* context */
{
  unsigned char bits[8];
  unsigned int index, padLen;

  /* Save number of bits */
  Encode (bits, context->count, 8);
  /* Pad out to 56 mod 64.
*/
  index = (unsigned int)((context->count[0] >> 3) & 0x3f);
  padLen = (index < 56) ? (56 - index) : (120 - index);
  MD5Update (context, PADDING, padLen);

  /* Append length (before padding) */
  MD5Update (context, bits, 8);
  /* Store state in digest */
  Encode (digest, context->state, 16);

  /* Zeroize sensitive information.
*/
  MD5_memset ((POINTER)context, 0, sizeof (*context));
}

/* MD5 basic transformation. Transforms state based on block.
 */
static void MD5Transform (state, block)
UINT4 state[4];
unsigned char block[64];
{
  UINT4 a = state[0], b = state[1], c = state[2], d = state[3], x[16];

  Decode (x, block, 64);

  /* Round 1 */
  FF (a, b, c, d, x[ 0], S11, 0xd76aa478); /* 1 */
  FF (d, a, b, c, x[ 1], S12, 0xe8c7b756); /* 2 */
  FF (c, d, a, b, x[ 2], S13, 0x242070db); /* 3 */
  FF (b, c, d, a, x[ 3], S14, 0xc1bdceee); /* 4 */
  FF (a, b, c, d, x[ 4], S11, 0xf57c0faf); /* 5 */
  FF (d, a, b, c, x[ 5], S12, 0x4787c62a); /* 6 */
  FF (c, d, a, b, x[ 6], S13, 0xa8304613); /* 7 */
  FF (b, c, d, a, x[ 7], S14, 0xfd469501); /* 8 */
  FF (a, b, c, d, x[ 8], S11, 0x698098d8); /* 9 */
  FF (d, a, b, c, x[ 9], S12, 0x8b44f7af); /* 10 */
  FF (c, d, a, b, x[10], S13, 0xffff5bb1); /* 11 */
  FF (b, c, d, a, x[11], S14, 0x895cd7be); /* 12 */
  FF (a, b, c, d, x[12], S11, 0x6b901122); /* 13 */
  FF (d, a, b, c, x[13], S12, 0xfd987193); /* 14 */
  FF (c, d, a, b, x[14], S13, 0xa679438e); /* 15 */
  FF (b, c, d, a, x[15], S14, 0x49b40821); /* 16 */

 /* Round 2 */
  GG (a, b, c, d, x[ 1], S21, 0xf61e2562); /* 17 */
  GG (d, a, b, c, x[ 6], S22, 0xc040b340); /* 18 */
  GG (c, d, a, b, x[11], S23, 0x265e5a51); /* 19 */
  GG (b, c, d, a, x[ 0], S24, 0xe9b6c7aa); /* 20 */
  GG (a, b, c, d, x[ 5], S21, 0xd62f105d); /* 21 */
  GG (d, a, b, c, x[10], S22,  0x2441453); /* 22 */
  GG (c, d, a, b, x[15], S23, 0xd8a1e681); /* 23 */
  GG (b, c, d, a, x[ 4], S24, 0xe7d3fbc8); /* 24 */
  GG (a, b, c, d, x[ 9], S21, 0x21e1cde6); /* 25 */
  GG (d, a, b, c, x[14], S22, 0xc33707d6); /* 26 */
  GG (c, d, a, b, x[ 3], S23, 0xf4d50d87); /* 27 */
  GG (b, c, d, a, x[ 8], S24, 0x455a14ed); /* 28 */
  GG (a, b, c, d, x[13], S21, 0xa9e3e905); /* 29 */
  GG (d, a, b, c, x[ 2], S22, 0xfcefa3f8); /* 30 */
  GG (c, d, a, b, x[ 7], S23, 0x676f02d9); /* 31 */
  GG (b, c, d, a, x[12], S24, 0x8d2a4c8a); /* 32 */

  /* Round 3 */
  HH (a, b, c, d, x[ 5], S31, 0xfffa3942); /* 33 */
  HH (d, a, b, c, x[ 8], S32, 0x8771f681); /* 34 */
  HH (c, d, a, b, x[11], S33, 0x6d9d6122); /* 35 */
  HH (b, c, d, a, x[14], S34, 0xfde5380c); /* 36 */
  HH (a, b, c, d, x[ 1], S31, 0xa4beea44); /* 37 */
  HH (d, a, b, c, x[ 4], S32, 0x4bdecfa9); /* 38 */
  HH (c, d, a, b, x[ 7], S33, 0xf6bb4b60); /* 39 */
  HH (b, c, d, a, x[10], S34, 0xbebfbc70); /* 40 */
  HH (a, b, c, d, x[13], S31, 0x289b7ec6); /* 41 */
  HH (d, a, b, c, x[ 0], S32, 0xeaa127fa); /* 42 */
  HH (c, d, a, b, x[ 3], S33, 0xd4ef3085); /* 43 */
  HH (b, c, d, a, x[ 6], S34,  0x4881d05); /* 44 */
  HH (a, b, c, d, x[ 9], S31, 0xd9d4d039); /* 45 */
  HH (d, a, b, c, x[12], S32, 0xe6db99e5); /* 46 */
  HH (c, d, a, b, x[15], S33, 0x1fa27cf8); /* 47 */
  HH (b, c, d, a, x[ 2], S34, 0xc4ac5665); /* 48 */

  /* Round 4 */
  II (a, b, c, d, x[ 0], S41, 0xf4292244); /* 49 */
  II (d, a, b, c, x[ 7], S42, 0x432aff97); /* 50 */
  II (c, d, a, b, x[14], S43, 0xab9423a7); /* 51 */
  II (b, c, d, a, x[ 5], S44, 0xfc93a039); /* 52 */
  II (a, b, c, d, x[12], S41, 0x655b59c3); /* 53 */
  II (d, a, b, c, x[ 3], S42, 0x8f0ccc92); /* 54 */
  II (c, d, a, b, x[10], S43, 0xffeff47d); /* 55 */
  II (b, c, d, a, x[ 1], S44, 0x85845dd1); /* 56 */
  II (a, b, c, d, x[ 8], S41, 0x6fa87e4f); /* 57 */
  II (d, a, b, c, x[15], S42, 0xfe2ce6e0); /* 58 */
  II (c, d, a, b, x[ 6], S43, 0xa3014314); /* 59 */
  II (b, c, d, a, x[13], S44, 0x4e0811a1); /* 60 */
  II (a, b, c, d, x[ 4], S41, 0xf7537e82); /* 61 */
  II (d, a, b, c, x[11], S42, 0xbd3af235); /* 62 */
  II (c, d, a, b, x[ 2], S43, 0x2ad7d2bb); /* 63 */
  II (b, c, d, a, x[ 9], S44, 0xeb86d391); /* 64 */

  state[0] += a;
  state[1] += b;
  state[2] += c;
  state[3] += d;

  /* Zeroize sensitive information.
*/
  MD5_memset ((POINTER)x, 0, sizeof (x));
}

/* Encodes input (UINT4) into output (unsigned char). Assumes len is
  a multiple of 4.
 */
static void Encode (output, input, len)
unsigned char *output;
UINT4 *input;
unsigned int len;
{
  unsigned int i, j;

  for (i = 0, j = 0; j < len; i++, j += 4) {
 output[j] = (unsigned char)(input[i] & 0xff);
 output[j+1] = (unsigned char)((input[i] >> 8) & 0xff);
 output[j+2] = (unsigned char)((input[i] >> 16) & 0xff);
 output[j+3] = (unsigned char)((input[i] >> 24) & 0xff);
  }
}

/* Decodes input (unsigned char) into output (UINT4). Assumes len is
  a multiple of 4.
 */
static void Decode (output, input, len)
UINT4 *output;
unsigned char *input;
unsigned int len;
{
  unsigned int i, j;

  for (i = 0, j = 0; j < len; i++, j += 4)
 output[i] = ((UINT4)input[j]) | (((UINT4)input[j+1]) << 8) |
   (((UINT4)input[j+2]) << 16) | (((UINT4)input[j+3]) << 24);
}

/* Note: Replace "for loop" with standard memcpy if possible.
 */

static void MD5_memcpy (output, input, len)
POINTER output;
POINTER input;
unsigned int len;
{
  unsigned int i;

  for (i = 0; i < len; i++)
 output[i] = input[i];
}

/* Note: Replace "for loop" with standard memset if possible.
 */
static void MD5_memset (output, value, len)
POINTER output;
int value;
unsigned int len;
{
  unsigned int i;

  for (i = 0; i < len; i++)
 ((char *)output)[i] = (char)value;
}

  /* Encode a 16-character MD5 digest in base-64 (25 characters) */

int cbf_md5digest_to64 (char *encoded_digest, const unsigned char *digest)
{
  static char basis_64 [] =

       "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

  int todo;
  
  if (!encoded_digest || !digest)
  
    return 0;
  

    /* Encode the 16 characters in base 64 */
    
  for (todo = 0; todo < 18; todo += 3)
  {
    encoded_digest [0] = basis_64 [((digest [todo + 0] >> 2) & 0x03f)];

    if (todo < 15)
    {
      encoded_digest [1] = basis_64 [((digest [todo + 0] << 4) & 0x030) |
                                     ((digest [todo + 1] >> 4) & 0x00f)];
      encoded_digest [2] = basis_64 [((digest [todo + 1] << 2) & 0x03c) |
                                     ((digest [todo + 2] >> 6) & 0x003)];
      encoded_digest [3] = basis_64 [((digest [todo + 2])      & 0x03f)];
    }
    else
    {
      encoded_digest [1] = basis_64 [((digest [todo + 0] << 4) & 0x030)];

      encoded_digest [2] = encoded_digest [3] = '=';
    }

    encoded_digest += 4;
  } 
  
  *encoded_digest  = '\0';

  return 0;
}

double CriterionPhi(double x, double A, double B, double C, double a, double b , double c)
{
  double Line;
  Line = ((A*A*b + C*(b*C - B*c) + A*B*(x - a))/(A*A + C*C));
  return Line;
}

double Polarization(double theta , double phi)
{
  double CorrIntent; 
  CorrIntent = 1.0*(1.0 - pow((sin(theta)*cos(phi)) , 2.0)) +
               0.0*(1.0 - pow((sin(phi)*sin(theta)) , 2.0));
  return CorrIntent;
}


double ThicknessCorrection(double theta , double phi, double thickness, double AbsoCoef, double psi)
{     
  double CorrIntent1, CorrIntent2, v , w, t, mu;
  
  t = thickness;
  mu = AbsoCoef;
  if (theta > 0)
    {
      v = sin(theta) * sin(phi) * cos(psi) + cos(theta)*sin(psi);
      w = - mu * t * ((v - sin(psi))/(v * sin(psi)));
      CorrIntent1 = ((exp(- mu * t / v))/w)*(exp(w) - 1.0);      
    }
  else CorrIntent1 = 1.0;
  if (CorrIntent1 > 0)
    {
      CorrIntent2 = CorrIntent1;
    }
  else CorrIntent2 = 1.0;
  return CorrIntent2;
}

double AirAbsorption(double r, double PixelSize, double AirAbso)
{
  double CorrIntent;
  CorrIntent = exp(-AirAbso*(r*PixelSize));
  return CorrIntent;
}

double DotProduct (double r , double theta , double phi, 
		       double detx, double dety, double detz)
{
  double CorrIntent, xpr , ypr , zpr , dotproduct , normal , l , pcnt, normdot,factor;
  xpr = r*sin(theta)*cos(phi);
  ypr = r*sin(theta)*sin(phi);
  zpr = r*cos(theta);
  normal = (sqrt(xpr*xpr + ypr*ypr + zpr*zpr));
  xpr = xpr/normal;
  ypr = ypr/normal;
  zpr = zpr/normal;
  dotproduct = (xpr*detx + ypr*dety + zpr*detz);
  if(dotproduct < 0)
    dotproduct = -dotproduct;
  /*
   l = (0.0081 / dotproduct);
   pcnt = exp(-20.0*l);
   factor = 0.2;
   pcnt = (1-factor)/(1-exp(log(factor)/dotproduct));
*/
  pcnt = dotproduct;
  CorrIntent = pcnt;
  return CorrIntent;
}

double PixelThickness (double r , double theta , double phi, 
		       double detx, double dety, double detz, double c)
{
  double CorrIntent, xpr , ypr , zpr , dotproduct , normal , l , pcnt, normdot,factor1,factor2,alpha;
  xpr = r*sin(theta)*cos(phi);   
  ypr = r*sin(theta)*sin(phi);
  zpr = r*cos(theta);
  normal = (sqrt(xpr*xpr + ypr*ypr + zpr*zpr));
  xpr = xpr/normal;
  ypr =  ypr/normal;
  zpr = zpr/normal;
  dotproduct = (xpr*detx + ypr*dety + zpr*detz);
  if(dotproduct < 0)
    dotproduct = -dotproduct;
  alpha = acos(dotproduct);
 

   
  /*14keV*/ pcnt = -0.28045958186777*alpha*alpha*alpha*alpha*alpha*alpha +
	      0.29342945094737*alpha*alpha*alpha*alpha +
	      0.17697385583177*alpha*alpha + 1.0;


  CorrIntent = pcnt;
  return CorrIntent;
}
