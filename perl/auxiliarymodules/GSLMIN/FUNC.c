#include <stdio.h>
#include <gsl/gsl_multimin.h>

static SV* ext_funname1;
static SV* ext_funname2;

static pdl* px;
static pdl* pgrad;
static SV* pxsv;
static SV* pgradsv;

static int ene;

double FF(int* n, double* x);
double DFF(int* n, double* x, double* grad);


double my_f (const gsl_vector *v, void *params);
void my_df (const gsl_vector *v, void *params, gsl_vector *df);
void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df);


static void my_delete_magic(pdl *p,int pa){
  p->data = 0;
}
 
void dump_pdl_kmk(pdl*p, FILE*f)
{
  fprintf(f,"magic: %lu state: %d, datatype: %d \n ",p->magicno,p->state,p->datatype);
}

pdl* haz_pdl(){
  pdl* t;
 
  int ndims;
  PDL_Indx *nelem;
  PDL_Indx *pdims;

  t = PDL->pdlnew();
  ndims = 1;
  pdims = (PDL_Indx *)  PDL->smalloc( (ndims) * sizeof(*pdims) );
  pdims[0] = (PDL_Indx) ene;
  PDL->setdims (t,pdims,ndims);
  t->datatype = PDL_D;
  PDL->allocdata (t);
  PDL->make_physical (t);
  t->state |= PDL_DONTTOUCHDATA | PDL_ALLOCATED;
  PDL->add_deletedata_magic(t,my_delete_magic,0);

  SV *ref;
  HV *stash = gv_stashpv("PDL",TRUE);
  SV *psv = newSViv(PTR2IV(t));
  t->sv = psv;

  return t;
}
 
//SV* haz_pdlsv(){
//  pdl* t;
//  SV* newref;
//
//  int ndims;
//  PDL_Indx *nelem;
//  PDL_Indx *pdims;
//
//  t = PDL->create(PDL_PERM);
//  ndims = 1;
//  pdims = (PDL_Indx *)  PDL->smalloc( (ndims) * sizeof(*pdims) );
//  pdims[0] = (PDL_Indx) ene;
//  PDL->setdims (t,pdims,ndims);
//  t->datatype = PDL_D;
//  PDL->allocdata (t);
//  PDL->make_physical (t);
//  
//  
//  SV *ref;
//  HV *stash = gv_stashpv("PDL",TRUE);
//  SV *psv = newSViv(PTR2IV(t));
//  t->sv = psv;
//  newref = newRV_noinc(t->sv);
//  (void)sv_bless(newref,stash);
//  return newref;
//}
 

double FF(int* n, double* x){

  double* res;
  double retval;
  int count;
  I32 ax ; 


  res = (double *) malloc(sizeof(double));
  px = haz_pdl();

  px->data = (void *) x;

  dSP;
  SV* funname;

  /* get function name on the perl side */
  funname = ext_funname1;

  ENTER;
  SAVETMPS;

  HV *stash = gv_stashpv("PDL",TRUE);
  pxsv = newRV_noinc(px->sv);
  (void) sv_bless(pxsv,stash);


  PUSHMARK(SP);

  XPUSHs(pxsv);

  PUTBACK;

  count=call_sv(funname,G_SCALAR);

  SPAGAIN; 
  SP -= count ;
  ax = (SP - PL_stack_base) + 1 ;

  if (count!=1)
    croak("error calling perl function\n");

  /* recover output value */
  *res = SvNV(ST(0));


  PUTBACK;
  FREETMPS;
  LEAVE;
  
  retval = *res;
  free(res);

  //  printf("FF: retval  %g \n", retval);

  return retval;
}



double DFF(int* n, double* xval, double* grad){

  double* res; double retval;
  double* xpass; int i;
  int count;
  I32 ax ; 


  res = (double *) malloc(sizeof(double));
  //  xpass = (double *) malloc( (*n) * sizeof(double));

  px = haz_pdl();
  pgrad = haz_pdl();

  px->data = (void *) xval;
  pgrad->data = (void *) grad;

  dSP;
  SV* funname;

  /* get function name on the perl side */
  funname = ext_funname2;

  ENTER;
  SAVETMPS;

  HV *stash = gv_stashpv("PDL",TRUE);
  pxsv = newRV_noinc(px->sv);
  (void) sv_bless(pxsv,stash);

  pgradsv = newRV_noinc(pgrad->sv);
  (void) sv_bless(pgradsv,stash);

  PUSHMARK(SP);

  XPUSHs(pxsv);

  //  XPUSHs(pgradsv);

  PUTBACK;

  //  count=call_sv(funname,G_SCALAR);

  //  printf("in DFF, about to call funnname \n");

  count=call_sv(funname,G_ARRAY);

  // printf("in DFF, outoffunname \n");

  SPAGAIN; 
  SP -= count ;
  ax = (SP - PL_stack_base) + 1 ;

  if (count!=1)
    croak("error calling perl function\n");

  /* recover output value */

  //   printf("in DFF, recovering pgrad \n");

   pgrad = PDL->SvPDLV(ST(0));
   PDL->allocdata(pgrad);
   PDL->make_physical(pgrad);
   PDL->make_physdims(pgrad);

   //   dump_pdl_kmk(pgrad,stderr);

   //   fprintf(stderr, "Type %d\n", pgrad->datatype);
   //   xpass  = (double *) pgrad->data;

   xpass  =  pgrad->data;

   //   printf("in DFF, about to assign xpass \n");

   for(i=0;i<ene;i++) {
     grad[i] =  xpass[i];
     //    printf("in DFF, output deriv is %g \n",grad[i]);
   }

   //  printf("DFF: grad (0)  %g \n", grad[0]);

   //  *res = SvNV(ST(0));


  PUTBACK;
  FREETMPS;
  LEAVE;
  
  // retval = *res;
  free(res);
  //  printf("in DFF, leaving \n");

  return 0.;
}

//
//void my_f (double* xfree, int  *nelem, double* out)
//{
//  int i;
//  float *xpass;
//  double retval;
//
////  printf ("in my_f, input is \n");
////  for (i=0;i<*nelem; i++) {
////     printf ("%g", xfree[i]);
////     //     xpass[i] = (float) (xfree[i]);
////  }
////  printf ("nelem: %d ", *nelem);
// 
//  retval = FF( nelem, xfree);
//
//  //  printf(" FF output is %g \n ",retval);
//
//  *out =  retval;
//
//  //  free(xpass);
//
//  //  printf(" my_f output is %g \n ",*out);
//
//  return;
//}
//
//
//void my_df (double* xfree, int  *nelem, double* out)
//{
//  int i;
//  float *xpass;
//  double retval;
//
////  printf ("in my_f, input is \n");
////  for (i=0;i<*nelem; i++) {
////     printf ("%g", xfree[i]);
////     //     xpass[i] = (float) (xfree[i]);
////  }
////  printf ("nelem: %d ", *nelem);
// 
//  retval = DFF( nelem, xfree);
//
//  //  printf(" FF output is %g \n ",retval);
//
//  *out =  retval;
//
//  //  free(xpass);
//
//  //  printf(" my_f output is %g \n ",*out);
//
//  return;
//}


//######################################################################


double
my_f (const gsl_vector *v, void *params)
{
  double *dp = (double *) params;
  int* nelem; int i;
  double* xfree;
  double retval;

  //  printf("testing my_f \n");

  nelem = (int *) malloc(sizeof(int));
  
  *nelem = (int) dp[0];

  //  printf("testing my_f : nelem is %d  \n", *nelem);
  
  xfree = (double *) malloc(*nelem * sizeof(double));

  
  for (i=0;i< *nelem;i++) {

    xfree[i] = gsl_vector_get(v, i);
    //    printf("testing my_f : i %d xfree %g \n", i, xfree[i]);
  }
  
  //  printf("testing my_f : i %d xfree %25.8e \n", 0, xfree[0]);

  retval = FF( nelem, xfree);
  //  printf(" FF output is %g \n ",retval);
  free(nelem);
  free(xfree);
  return retval ;
}

/* The gradient of f, df = (df/dx, df/dy). */
void
my_df (const gsl_vector *v, void *params,
       gsl_vector *df)
{
  double *dp = (double *)params;
  int* nelem; int i;
  double* xfree;
  double* grad;
  double retval;

  //  printf("testing my_df \n");
  
  nelem = (int *) malloc(sizeof(int));
  *nelem = (int) dp[0];

  //printf("testing my_df : nelem is %d  \n", *nelem);

  xfree = (double *) malloc(*nelem * sizeof(double));
  grad = (double *) malloc(*nelem * sizeof(double));
  
  
  for (i=0;i< *nelem;i++) {
    xfree[i] = gsl_vector_get(v, i);
    grad[i] = gsl_vector_get(v, i);
    // printf("testing my_df input vector: i  %d xfree %g \n", i, xfree[i]);
  }
  
  DFF(nelem, xfree, grad);
  
  for (i=0;i< *nelem;i++) {
    gsl_vector_set(df, i, grad[i]);
    // printf("testing my_df output deriv: i %d xfree %g \n", i, xfree[i]);
  }
  

  free(nelem);
  free(xfree);
  free(grad);
  return;
  
}

/* Compute both f and df together. */
void
my_fdf (const gsl_vector *x, void *params,
	double *f, gsl_vector *df)
{
  *f = my_f(x, params);
  my_df(x, params, df);
}


void conjgrad (double *xfree, int  *nelem, double *out, double *stepsize, double *linmintol, double *gradtol) 
{
  size_t iter = 0;
  int status; int i;
  double par[1] = { *nelem };
  double size;
  double fprev = +1E30; // ADDITION BY SIMON TO CONTROL CONVERGENCE BY VALUE RATHER THAN GRADIENT
  double reldiff = 0; // ADDITION BY SIMON TO CONTROL CONVERGENCE BY VALUE RATHER THAN GRADIENT

  //printf ("IN TEST: nelem = %d \n", *nelem);
//  for (i=0; i<*nelem; i++) {
//    printf ("IN TEST: %g \n", xfree[i]);
//  }
//

  printf ("conjgrad input params : stepsize %g linmintol %g gradtol %g \n", *stepsize, *linmintol, *gradtol); 
  
  

  const gsl_multimin_fdfminimizer_type *T;

//  printf ("ANCHOR \n");
//  printf ("IN TEST2: nelem = %d \n", *nelem);
//
  gsl_multimin_fdfminimizer *s;
  
//  printf ("ANCHOR \n");
//  printf ("IN TEST3: nelem = %d \n", *nelem);
//

  gsl_vector *x;

  //  printf ("ANCHOR \n");
  // printf ("IN TEST4: nelem = %d \n", *nelem);

  gsl_multimin_function_fdf my_func;
  
  //printf ("ANCHOR \n");
  // printf ("IN TEST5: nelem = %d \n", *nelem);

  x = gsl_vector_alloc (*nelem);

  // printf ("ANCHOR \n");
  // printf ("IN TEST6: nelem = %d \n", *nelem);

  my_func.f = &my_f;

  // printf ("ANCHOR \n");
  // printf ("IN TEST7: nelem = %d \n", *nelem);

  my_func.df = &my_df;
  my_func.fdf = &my_fdf;

  //  printf ("ANCHOR \n");
  //  printf ("IN TEST8: nelem = %d \n", *nelem);
  //
  //  my_func.n = *nelem;
  //my_func.params = *nelem;



  my_func.n = *nelem;



  my_func.params = par;
  
  //  printf ("nelem: %d \n", my_func.n);
  //  printf ("in conjgrad, input is \n");
  
  //  for (i=0; i< **nelem; i++) {
  //    //     printf ("%g", xfree[i]);
  //     printf ("%d ", i);
  //  }



  //printf ("IN TEST: nelem = %d \n", *nelem);
  //  for (i=0; i<*nelem; i++) {
  //   printf ("IN TEST: %g \n", xfree[i]);
  // }


  for (i=0;i<*nelem; i++) {
    gsl_vector_set (x, i, xfree[i]);
  }


  //  T = gsl_multimin_fdfminimizer_conjugate_pr;
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, *nelem);
  

  gsl_multimin_fdfminimizer_set (s, &my_func, x, *stepsize, *linmintol);
 
  // my_f (const gsl_vector *v, void *params)

	  
  do {



    iter++;

    printf( "GSL FDF ITERATE \n");
    printf ("conjgrad active  params : stepsize %g linmintol %g gradtol %g \n", *stepsize, *linmintol, *gradtol); 

    //    status = gsl_multimin_fdfminimizer_iterate (s);
    (s->type->iterate) (s->state, s->fdf, s->x, &(s->f), s->gradient, s->dx);

    printf("done iterate\n");
    // from fdfminimizer.c :
    // int
    // gsl_multimin_fdfminimizer_iterate (gsl_multimin_fdfminimizer * s)
    // {
    //   return (s->type->iterate) (s->state, s->fdf, s->x, &(s->f), s->gradient, s->dx);
    // }
    // 
    // 


    //    printf( "GSL FDF ITERATE DONE, status %d \n",status);
    
    
    //    printf( "TESTING GRADIENT TOL %g   \n", *gradtol);
    status = gsl_multimin_test_gradient (s->gradient, *gradtol);
    

    if (status == GSL_CONTINUE) {
      //      printf("TESTING VALUE FOR TOL %e \n",1E-6);
      if (fprev == 0) {
	if (s->f == 0) {
	  status = GSL_SUCCESS;  // ADITION BY SIMON TO CONTROL CONVERGENCE BY VALUE
	  printf("CONVERGENCE TO ZERO BY VALUE - not by gradient ! \n");
	}
      } else {
	reldiff = fabs( (s->f - fprev) / fprev);
	if (reldiff < 1E-8) { 
	  status = GSL_SUCCESS;
	  printf("CONVERGENCE BY VALUE - not by gradient, as the GSL default ! \n");
	}
      }
    }

    fprev = s->f;

    // -- Function: gsl_vector * gsl_multimin_fdfminimizer_gradient (const
    //     gsl_multimin_fdfminimizer * S)


    //    status = gsl_multimin_test_size (size, 1e-2);




    //    //    size = my_f(x, par);
    //    size = s->f;
    //    printf( "TESTING SIZE %g \n", size);
    //    status = gsl_multimin_test_size (size, *gradtol);

    
    //    printf( "GSLSTATUS  %d cont: %d success: %d  \n",status, GSL_CONTINUE, GSL_SUCCESS);

	
    //    if (status == GSL_SUCCESS)
    //  printf ("Minimum found at:\n");

    
    // for (i=0;i<*nelem; i++) {
    //  printf (" %.5f ", gsl_vector_get (s->x, i));
    // }

    printf ("iter %5d size %20.15e\n", iter, s->f); // s->f gives current function value

//     for (i=0;i<*nelem; i++) {
//       printf("gradient: %d %g  ",i, gsl_vector_get(s->gradient, i));
//     }

    for (i=0;i<*nelem; i++) {
      xfree[i] = gsl_vector_get (s->x, i);
    }
    
  } while (status == GSL_CONTINUE && iter < 50);
  //  while (iter < 1000000);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
  
  return ;
}
