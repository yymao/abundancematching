#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define GAUSS_CACHE_SIZE 4097
#define GAUSS_CACHE_CENTER 2048 /* (4097-1)/2 */
#define GAUSS_CACHE_MAX 8
#define GAUSS_CACHE_INV_STEP 256 /* 2048/8 */

double gaussian_cache(double dm, double scatter, double * g_cache)
{
  int b;
  dm = dm/scatter*GAUSS_CACHE_INV_STEP + GAUSS_CACHE_CENTER;
  b = dm;
  if (b<0 || b>=GAUSS_CACHE_SIZE-1) return 0;
  dm -= b;
  return (g_cache[b] + dm*(g_cache[b+1]-g_cache[b]))/scatter;
}

void convolved_fit(double * af_key, double * af_val, int num_af_points, 
    double * smm, double * mf, int MASS_BINS, double scatter, 
    int repeat, double sm_step)
{
  if (!repeat || !scatter) return;

  int i, j, k, jstart, jend;
  double sm, fid_sm; 
  double nh, new_nh, lnh, lnnh, mnh, sum;

  double * new_smm = (double *) malloc(MASS_BINS*sizeof(double));
  double * smm_conv = (double *) malloc(MASS_BINS*sizeof(double));
  double * g_cache = (double *) malloc(GAUSS_CACHE_SIZE*sizeof(double));

  for (i=0; i<GAUSS_CACHE_SIZE; i++) {
    sm = ((double)i)/GAUSS_CACHE_INV_STEP - GAUSS_CACHE_MAX;
    g_cache[i] = exp(-0.5*sm*sm)/sqrt(2.0*3.1415926535);
  }

  for (k=0; k<repeat; ++k) {
    jstart = 0;
    jend = MASS_BINS;
    fid_sm = af_key[num_af_points-1]+1;
    nh = pow(10, af_val[num_af_points-1])
        * (af_key[num_af_points-1] - af_key[num_af_points-2]);

    new_smm[0] = 0;
    for (i=1; i<MASS_BINS; i++) {
      new_smm[i] = 0;
      if (smm[i] && !smm[i-1]) jstart = i;
      if (!smm[i] && smm[i-1]) jend = i;
    }
    
    mnh = 0;
    while (mnh<nh) {
      sum = 0;
      for (j=jstart; j<jend; j++) {
        smm_conv[j] = gaussian_cache(fid_sm-smm[j], scatter, g_cache);
        sum += smm_conv[j]*mf[j];
      }
      mnh += sum * sm_step;
      if (!isfinite(mnh)) mnh = 0;
      for (j=jstart; j<jend; j++) {
        new_smm[j] += smm_conv[j]*sm_step*af_key[num_af_points-1]; 
      }
      fid_sm -= sm_step;
    }
    
    lnh = log10(nh);
    for (i=num_af_points - 1; i>0; i--) {
      new_nh = nh + (pow(10, af_val[i-1])-pow(10, af_val[i]))
          * (af_key[i] - af_key[i-1]) 
          / (log(10)*(af_val[i-1]-af_val[i]));
      lnnh = log10(new_nh);
      while (mnh<new_nh && fid_sm>=0) {
        sum = 0;
        for (j=jstart; j<jend; j++) {
          smm_conv[j] = gaussian_cache(fid_sm-smm[j], scatter, g_cache);
          sum += smm_conv[j]*mf[j];
        }
        mnh += sum * sm_step;
        sm = af_key[i-1] + (log10(mnh)-lnnh)/(lnh-lnnh) 
            * (af_key[i]-af_key[i-1]);
        for (j=jstart; j<jend; j++) new_smm[j] += smm_conv[j]*sm*sm_step;
        fid_sm -= sm_step;
      }
      nh = new_nh;
      lnh = lnnh;
    }

    for (j=jstart; j<jend; j++) smm[j] = new_smm[j];
  }

  free(new_smm);
  free(smm_conv);
  free(g_cache);
}

