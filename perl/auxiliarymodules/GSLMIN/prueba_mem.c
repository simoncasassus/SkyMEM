#include <stdio.h>
#include <gsl/gsl_multimin.h>

/* Paraboloid centered on (dp[0],dp[1]) */

double
my_f (const gsl_vector *v, void *params)
{
  double x, y;
  double *dp = (double *) params;
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  printf("f\n");
  return (x - dp[0]) * (x - dp[0]) +
    (y - dp[1]) * (y - dp[1]) + 3.0;
}

/* The gradient of f, df = (df/dx, df/dy). */
void
my_df (const gsl_vector *v, void *params,
       gsl_vector *df)
{
  double x, y;
  double *dp = (double *)params;
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
  
  gsl_vector_set(df, 0, 2.0 * (x - dp[0]));
  gsl_vector_set(df, 1, 2.0 * (y - dp[1]));
}

/* Compute both f and df together. */
void
my_fdf (const gsl_vector *x, void *params,
	double *f, gsl_vector *df)
{
  *f = my_f(x, params);
  my_df(x, params, df);
}


int main (void)
{
  size_t iter = 0;
  int status;
  /***********************************************************/

  long *seed = 11;
  int i, k;
  int N = 128;
  double im[N][N];
  double free[N*N];
  double signoise = 0.1;
  double sigwidth = 5;

  printf("INPUT *seed %d \n", (int) *seed);

  for (i = 0; i < N; i++) {      /* Loop through samples */
    for (k = 0; k < N; k++) {    /* Loop through channels */
      
      im[i][k] = exp(-0.5*(SQR(i-N/2)+SQR(j-N/2))/SQR(sigwidth)) +
	exp(-0.5*(SQR(i-7-N/2)+SQR(j-7-N/2))/SQR(sigwidth));
      
      

      im[i][k]  += gasdev(seed)*signoise;   // DEVTEST
      }    /* weights */
    }
  }

  /***********************************************************/



  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  
  /* Position of the minimum (1,2). */
  double par[2] = { 2.0, 0.0 };
  
  gsl_vector *x;
  gsl_multimin_function_fdf my_func;
  
  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n = 2;
  my_func.params = &par;


  /* Starting point, x = (5,7) */
  
  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 12.2);
  gsl_vector_set (x, 1, 20000.0);
  
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, 2);
  
  gsl_multimin_fdfminimizer_set (s, &my_func, x, 1e-6, 1e-4);
  
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status)
      break;
    
    status = gsl_multimin_test_gradient (s->gradient, 1e-4);
	
    if (status == GSL_SUCCESS)
      printf ("Minimum found at:\n");
    
    printf ("%5d %.5f %.5f %10.5f\n", iter,
	    gsl_vector_get (s->x, 0),
	    gsl_vector_get (s->x, 1),
	    s->f);
    
  }
  while (status == GSL_CONTINUE && iter < 1007);
  
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
  
  return 0;
}

