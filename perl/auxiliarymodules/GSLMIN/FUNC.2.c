static SV* ext_funname1;

static pdl* px;
static SV* pxsv;
static int ene;

float FF(int* n, float* x);

static void my_delete_magic(pdl *p,int pa){
  p->data = 0;
}
 
pdl* haz_pdl(){
  pdl* t;
 
  int ndims;
  PDL_Long *nelem;
  PDL_Long *pdims;

  t = PDL->pdlnew();
  ndims = 1;
  pdims = (PDL_Long *)  PDL->smalloc( (ndims) * sizeof(*pdims) );
  pdims[0] = (PDL_Long) ene;
  PDL->setdims (t,pdims,ndims);
  t->datatype = PDL_F;
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
  

float FF(int* n, float* x){

  float* res;
  int count;
  I32 ax ; 

  res = (float *) malloc(sizeof(float));
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
  
  return *res;
}


void my_f (double* xfree, int  *nelem, double* out)
{
  int i;
  float *xpass;
  //  (*out) = FF( (*nelem), (*xfree));

  xpass = (float*) malloc ( (*nelem) * sizeof(float)  ); 

  printf ("in my_f, input is \n");
  for (i=0;i<*nelem; i++) {
     printf ("%g", xfree[i]);
     xpass[i] = (float) (xfree[i]);
  }
  printf ("nelem: %d ", *nelem);
 
  (*out) =  FF( nelem, xpass);
  return;
}

