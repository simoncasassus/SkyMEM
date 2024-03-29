/*
  
  cc -g -o a.out  wcs_call.c libwcs/libwcs.a -lm

#include <stdio.h>
#include <string.h>
 
int main()
{
   char s[13] = "Hola a todos";
 
   printf( "s=%s\n", s );
   printf( "strlen(s) = %d\n", strlen( s ) );
 
   return 0;
}

    * length EXPR

    * length

      Returns the length in characters of the value of EXPR. If EXPR is omitted, returns length of $_ . Note that this cannot be used on an entire array or hash to find out how many elements these have. For that, use scalar @array and scalar keys %hash respectively.

      Note the characters: if the EXPR is in Unicode, you will get the number of characters, not the number of bytes. To get the length in bytes, use do { use bytes; length(EXPR) } , see bytes.


 */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "libwcs/wcs.h"
#include "libwcs/fitsfile.h"

extern struct WorldCoor *GetFITSWCS ();	/* Read WCS from FITS or IRAF header */
extern char *GetFITShead();
extern struct WorldCoor *GetWCSFITS ();	/* Read WCS from FITS or IRAF file */





int
proj (int direction, char* fn, int dohead, int n, int m, double imra[n][m], double imdec[n][m], double imx[n][m], double imy[n][m])
{
  int verbose = 0;		/* verbose/debugging flag */
  
  double x, y, ra, dec; // ra0, dec0;
  
  struct WorldCoor *wcs;
   
  char *header;

  //  char *header_pass = calloc(80*500,sizeof(char)); 

  char *fn_pass = calloc(50,sizeof(char)); 


  double cra, cdec, dra, ddec, secpix;
  int wp, hp;
  int sysout = 0;
  double eqout = 0.0;
  int offscale;
  
  int i =0;
  int j =0;

  


  if (!dohead) {
    header = GetFITShead (fn, verbose);

    //    printf("FILE PASS in wcs_call.c >>>>>>>>>>>>>>>>>>>\n!!!!%s!!!!\n",header);


    wcs = GetFITSWCS (fn,header,verbose,&cra,&cdec,&dra,&ddec,&secpix,
		    &wp, &hp, &sysout, &eqout);

  } else {


    strcpy(fn_pass,"dum.fits");



    //    char *fntest = calloc(12,sizeof(char)); 
    // strcpy(fntest,"ngc6302.fits");



    char *fntest = calloc(50,sizeof(char)); 

    //    strcpy(fntest,"nvss_polaris.fits");
    strcpy(fntest,"/home/simon/common/perl/pdlC/WCS/ngc6302.fits");


    
    //    char * fntest = "ngc6302.fits";

    header = GetFITShead (fntest, verbose);

    //    printf("HEADER PASS in wcs_call.c >>>>>>>>>>>>>>>>>>>\n%s\n",fn);
    //    printf("^^^^^^^^^HEADER PASS in wcs_call.c ^^^^^^^^^^^^^^^^^^^^\n");
//     printf("REAL HEAD  in wcs_call.c >>>>>>>>>>>>>>>>>>>\n!!!!%s!!!!\n",header);
//     printf("^^^^^^^^^REAL HEADER  in wcs_call.c ^^^^^^^^^^^^^^^^^^^^\n");
    //    wcs = GetFITSWCS (fn_pass,fn,verbose,&cra,&cdec,&dra,&ddec,&secpix,
    //		      &wp, &hp, &sysout, &eqout);


    wcs = GetFITSWCS (fn_pass,fn,verbose,&cra,&cdec,&dra,&ddec,&secpix,
		      &wp, &hp, &sysout, &eqout);

    //    printf("HEADER PASS PASSED\n");

   free(fn_pass);
   free(fntest);


  }



  if (direction == -1) {

    for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++) {

	
	ra = imra[i][j];
	dec = imdec[i][j];

	wcs2pix (wcs, ra, dec,  &x, &y, &offscale);
	
	imx[i][j] = x;
	imy[i][j] = y;

	
      }
    }
    
  } else if (direction == 1) {

    
    for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++) {
	
	x = imx[i][j];
	y = imy[i][j];
    
	
	pix2wcs (wcs, x, y,  &ra, &dec);
	
	
     	imra[i][j] = ra;
	imdec[i][j] = dec;
		
      }
    }

  } else {

    printf ("abs(direction) != 1 !\n");

    exit;
  }



  return 0;
}

//  gcc -o ttest wcs_call.c -L./libwcs/ -lwcs -lm

int main (void) 
{
  
  int direction = -1;

  int fnlength = 17;
  char *fn = calloc(fnlength, sizeof(char)); 


  int hdlength = (100*80);
  char *hd = calloc(hdlength, sizeof(char)); 


  int dohead = 0;

  int n=1;
  int m=1; 

  int i;

  double imra[n][m];
  double imdec[n][m];
  double imx[n][m];
  double imy[n][m];
  

  imra[0][0] =  167.481458333333 ;
  imdec[0][0] =  88.94575 ;


  strcpy( fn,"nvss_polaris.fits"); 
  //strcpy( fn,""); 
  
  proj (direction, fn,  dohead,  n,  m, imra, imdec, imx, imy);
  
  //  printf("imx %g\n",imx[0][0]);
  //printf("imy %g\n",imy[0][0]);

      
  strcpy (hd,"SIMPLE  =                     1 /                                               BITPIX  =                   -32 /                                               NAXIS   =                     2 /                                               NAXIS1  =                   300 /                                               NAXIS2  =                   300 /                                               CRVAL1  =    37.954514999999994 /                                               CRVAL2  =     89.26410900000005 /                                               RADESYS =                   FK5 /                                               EQUINOX =                2000.0 /                                               CTYPE1  =              RA---TAN /                                               CTYPE2  =              DEC--TAN /                                               CRPIX1  =                 150.5 /                                               CRPIX2  =                 150.5 /                                               CDELT1  =     -0.01666666666667 /                                               CDELT2  =  0.016666666666666666 /                                               COMMENT =                       /                                               SURVEY  =                  NVSS /                                               TELESCOP=                   VLA /                                               INSTRUME=                L-BAND /                                               OBSERVER=              NVSS GRP /                                               BUNIT   =               JY/BEAM /                                               HISTORY                                                                         END                                                                             ");


  
  char *header;
  int verbose = 0;		/* verbose/debugging flag */

  dohead = 1;

  proj (direction, hd,  dohead,  n,  m, imra, imdec, imx, imy);
  
  //  printf("imx %g\n",imx[0][0]);
  //printf("imy %g\n",imy[0][0]);
  
  free(hd);

  return 0;


}
