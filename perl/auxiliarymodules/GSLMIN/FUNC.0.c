static SV* ext_funname1;
static SV* ext_funname2;
//static gsl_multimim_function_fdf F;


static pdl* px;
static SV* pxsv;
static int ene;

double FF(int* n, double* x);
void DFF(int* n, double* x, double* df);

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
  

double FF(int* n, double* x){

  double* res;
  int count;
  I32 ax ; 

  res = (double*) malloc(sizeof(float));
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


void  DFF(int* n, double* x, double* df){

  double* res;
  int count;
  I32 ax ; 

  res = (double*) malloc(sizeof(float));
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






void my_f (double* xfree, int nelem, double* out)
{
  (*out) = FF(nelem,xfree);
  return;
}

/* The gradient of f, df = (df/dx, df/dy). */
void
my_df (double* xfree, int nelem,
       double* deriv)
{


  DFF (nelem, xfree, deriv);


  return;

}
