#include <R.h>
#include <Rmath.h>
#include <math.h>

double phi(double x) {
  double out;
  if ((x >= 0) && (x < 1)) {
    out = 1.0;
  } else {
    out = 0.0;
  }
  return(out);
}

double psi(double x) {
  double out;
  if ((x >= 0) && (x < 0.5)) {
    out = 1.0;
  } else {
    if ((x >= 0.5) && (x < 1.0)) {
      out = -1.0;
    } else {
      out = 0.0;
    }
  }
  return(out);
}

double thresholding(double x, double u) {
  double out;
  if (x > 0) {
    out = x - u;
    if (out < 0) {
      out = 0.0;
    }
  } else {
    if (x < 0) {
      out = x + u;
      if (out > 0) {
	out = 0.0;
      }
    } else {
      out = 0.0;
    }
  }
  return(out);
}

double kmax(double *x, int k) {
  int i;
  double out;

  out = x[0];
  for (i = 0; i < k; i++) {
    if (out < x[i]) {
      out = x[i];
    }
  }
  return(out);
}

int maxloc(int n, double *x) {
  int i, max_loc;
  double xmax;

  xmax = x[0];
  max_loc = 0;
  for (i = 0; i < n; i++) {
    if (xmax < x[i]) {
      xmax = x[i];
      max_loc = i;
    }
  }
  return(max_loc);
}

int minloc(int n, double *x) {
  int i, min_loc;
  double xmin;

  xmin = x[0];
  min_loc = 0;
  for (i = 0; i < n; i++) {
    if (xmin > x[i]) {
      xmin = x[i];
      min_loc = i;
    }
  }
  return(min_loc);
}
 

void quicksort(double *numbers, int first, int last){
  int i, j, pivot;
  double temp;

  if (first < last){
    pivot = first;
    i = first;
    j = last;

    while(i<j){

      while(numbers[i] <= numbers[pivot] && i < last)
	i++;

      while(numbers[j] > numbers[pivot])
	j--;
      
      if(i < j){
	
	temp = numbers[i];
	numbers[i] = numbers[j];
	numbers[j] = temp;
      }
    }

    temp = numbers[pivot];
    numbers[pivot] = numbers[j];
    numbers[j] = temp;
    quicksort(numbers, first, j - 1);
    quicksort(numbers, j + 1, last);

  }
}

double positive(double x) {
  double out;
  if (x > 0) {
    out = x;
  } else {
    out = 0.0;
  }
  return(out);
}

double sign(double x) {
  double out;
  if (x > 0) {
    out = 1.0;
  } else {
    if (x < 0) {
      out = -1.0;
    } else {
      out = 0.0;
    }
  }
  return(out);
}

void wav1d_dm(int *n_in, int *J_in, double *dm) {
  int n = n_in[0];
  int J = J_in[0];
  int i, j, k;

  double *x;

  x = (double *)malloc(n * sizeof(double));

  for (i = 0; i < n; i++) {
    x[i] = (1.0 * i)/(1.0 * n);
    dm[i] = phi(x[i]);
  }

  for (j = 0; j < J; j++) {
    for (k = 0; k < pow(2, j); k++) {
      for (i = 0; i < n; i++) {
	dm[((int)(round(pow(2, j))) + k)*n + i] = pow(2, j/2.0) * psi(pow(2, j) * x[i] -  k);
      }
    }
  }

  free(x);
}

void wav2d_coef(int *n1_in, int *n2_in, double *Z, int *J1_in, int *J2_in, double *wcoef) {
  double *dm1, *dm2;
  int n1 = n1_in[0];
  int n2 = n2_in[0];
  int J1 = J1_in[0];
  int J2 = J2_in[0];
  int i1, i2, p1, p2, NC1, NC2;

  NC1 = (int)round(pow(2, J1));
  dm1 = (double *)malloc(n1 * NC1 * sizeof(double));
  wav1d_dm(n1_in, J1_in, dm1);

  NC2 = (int)round(pow(2, J2));
  dm2 = (double *)malloc(n2 * NC2 * sizeof(double));
  wav1d_dm(n2_in, J2_in, dm2);

  for (p1 = 0; p1 < NC1; p1++) {
    for (p2 = 0; p2 < NC2; p2++) {

      wcoef[p1 * NC2 + p2] = 0.0;

      for (i1 = 0; i1 < n1; i1++) {
	for (i2 = 0; i2 < n2; i2++) {

	  wcoef[p1 * NC2 + p2] = wcoef[p1 * NC2 + p2] + dm1[p1 * n1 + i1] * dm2[p2 * n2 + i2] * Z[i1 * n2 + i2];

	}
      }

      wcoef[p1 * NC2 + p2] = wcoef[p1 * NC2 + p2] /(1.0 * n1 * n2);
      
    }
  }

  free(dm1);
  free(dm2);

}

void wav2d_img(int *n1_in, int *n2_in, int *J1_in, int *J2_in, double *wcoef, double *thresh_in, double *img) {
  double *dm1, *dm2;
  int n1 = n1_in[0];
  int n2 = n2_in[0];
  int J1 = J1_in[0];
  int J2 = J2_in[0];
  int i1, i2, p1, p2, NC1, NC2;

  double thresh = thresh_in[0];

  NC1 = (int)(round(pow(2, J1)));
  dm1 = (double *)malloc(n1 * NC1 * sizeof(double));
  wav1d_dm(n1_in, J1_in, dm1);

  NC2 = (int)round(pow(2, J2));
  dm2 = (double *)malloc(n2 * NC2 * sizeof(double));
  wav1d_dm(n2_in, J2_in, dm2);

  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {

      img[i1 * n2 + i2] = 0.0;

      for (p1 = 0; p1 < NC1; p1++) {
      	for (p2 = 0; p2 < NC2; p2++) {

      	  img[i1 * n2 + i2] = img[i1 * n2 + i2] + dm1[p1 * n1 + i1] * dm2[p2 * n2 + i2] * thresholding(wcoef[p1 * NC2 + p2], thresh);
  
      	}
      }

    }
  }

  free(dm1);
  free(dm2);

}

void ALO_phase1(int *p_in, double *sig2_in, double *X, double *W) {
  int p = p_in[0];
  double sig2 = sig2_in[0];
  double *Xabs;
  int b, q, iq, shrink_indx;
  double num, den, shrink, temp;

  double const small_number = pow(10, -37);

  Xabs = (double *)malloc(p * sizeof(double));

  for (b = 0; b < p; b++) {
    Xabs[b] = fabs(X[b]);
  }

  quicksort(Xabs, 0, p - 1);

  if (Xabs[0] < small_number) {
    error("X contains a 0 value. \n");
  }

  for (b = 0; b < (p - 1); b++) {
    if ((Xabs[b + 1] - Xabs[b]) < small_number) {
      error("X contains ties. \n");
    }
  }

  for (q = 0; q < (p - 1); q++) {
    shrink_indx = (p - 1) - (q + 1);
    shrink = Xabs[shrink_indx];
    num = 0.0;
    den = 0.0;
    for (iq = (shrink_indx + 1); iq < p; iq++) {
      temp = Xabs[iq] * Xabs[iq] -shrink * shrink;
      num = num + temp;
      den = den + temp * temp / (Xabs[iq] * Xabs[iq]);
    }
    W[q] = (num * num / sig2) / den;
  }

  W[p - 1] = 0.0;
  for (b = 0; b < p; b++) {
    W[p - 1] = W[p - 1] + Xabs[b] * Xabs[b];
  }
  W[p - 1] = W[p - 1]/sig2;

  free(Xabs);
}

double SCAD_1dmu(double std, double x, double zeta, double a) {
  double out;
  double xabs, xsign;

  xabs = fabs(x);
  xsign = sign(x);
  if (xabs <= (2 * zeta * std)) {
    out = xsign * positive(xabs - std * zeta);
  } else {
    if (xabs <= (a * zeta * std)) {
      out = ((a - 1) * x - xsign * a * zeta * std)/(a - 2);
    } else {
      out = x;
    }
  }
  return(out);
}

void SCAD_phase1(int *p_in, double *sig2_in, double *a_in, double *X, double *V) {
  int p = p_in[0];
  double sig2 = sig2_in[0];
  double a = a_in[0];
  double *Xabs, *MU;
  int b, q, shrink_indx;
  double num, den, shrink;

  double const small_number = pow(10, -37);

  Xabs = (double *)malloc(p * sizeof(double));

  for (b = 0; b < p; b++) {
    Xabs[b] = fabs(X[b]);
  }

  quicksort(Xabs, 0, p - 1);

  if (Xabs[0] < small_number) {
    error("X contains a 0 value. \n");
  }

  for (b = 0; b < (p - 1); b++) {
    if ((Xabs[b + 1] - Xabs[b]) < small_number) {
      error("X contains ties. \n");
    }
  }

  /* Compute the SCAD charting staistics */

  MU = (double *)malloc(p * sizeof(double));

  for (q = 0; q < (p - 1); q++) {
    shrink_indx = (p - 1) - (q + 1);
    shrink = Xabs[shrink_indx];

    /* Compute the SCAD estimate */
    for (b = 0; b < p; b++) {
      MU[b] = SCAD_1dmu(1.0, X[b], shrink, a);
    }

    num = 0.0;
    den = 0.0;
    for (b = 0; b < p; b++) {
      num = num + X[b] * MU[b];
      den = den + MU[b] * MU[b];
    }
    V[q] = (num * num / sig2) / den;
  }

  V[p - 1] = 0.0;
  for (b = 0; b < p; b++) {
    V[p - 1] = V[p - 1] + X[b] * X[b];
  }
  V[p - 1] = V[p - 1] / sig2;

  free(Xabs);
  free(MU);
     
}

void ALO_ICstat(int *n_IC_in, int *p_in, double *sig2_in, double *IC_bar, double *IC_std) {
  int n_IC = n_IC_in[0];
  int p = p_in[0];
  int j, b;
  double *X, *W, *W_mat;

  X = (double *)malloc(p * sizeof(double));
  W = (double *)malloc(p * sizeof(double));
  W_mat = (double *)malloc(n_IC * p * sizeof(double));

  GetRNGstate();

  for (j = 0; j < n_IC; j++) {

    for (b = 0; b < p; b++) {
      X[b] = sqrt(sig2_in[0]) * rnorm(0.0, 1.0);
    }

    ALO_phase1(p_in, sig2_in, X, W);

    for (b = 0; b < p; b++) {
      W_mat[j * p + b] = W[b];
    }

  }

  PutRNGstate();

  free(X);
  free(W);

  for (b = 0; b < p; b++) {
    IC_bar[b] = 0.0;
    IC_std[b] = 0.0;
    for (j = 0; j < n_IC; j++) {
      IC_bar[b] = IC_bar[b] + W_mat[j * p + b];
      IC_std[b] = IC_std[b] + W_mat[j * p + b] * W_mat[j * p + b];
    }
    IC_bar[b] = IC_bar[b] / (1.0 * n_IC);
    IC_std[b] = sqrt(IC_std[b] / (1.0 * n_IC) - IC_bar[b] * IC_bar[b]);
  }

  free(W_mat);

}

void SCAD_ICstat(int *n_IC_in, int *p_in, double *sig2_in, double *a_in, double *IC_bar, double *IC_std) {
  int n_IC = n_IC_in[0];
  int p = p_in[0];
  int j, b;
  double *X, *V, *V_mat;

  X = (double *)malloc(p * sizeof(double));
  V = (double *)malloc(p * sizeof(double));
  V_mat = (double *)malloc(n_IC * p * sizeof(double));

  GetRNGstate();

  for (j = 0; j < n_IC; j++) {

    for (b = 0; b < p; b++) {
      X[b] = sqrt(sig2_in[0]) * rnorm(0.0, 1.0);
    }

    SCAD_phase1(p_in, sig2_in, a_in, X, V);

    for (b = 0; b < p; b++) {
      V_mat[j * p + b] = V[b];
    }

  }

  PutRNGstate();

  free(X);
  free(V);

  for (b = 0; b < p; b++) {
    IC_bar[b] = 0.0;
    IC_std[b] = 0.0;
    for (j = 0; j < n_IC; j++) {
      IC_bar[b] = IC_bar[b] + V_mat[j * p + b];
      IC_std[b] = IC_std[b] + V_mat[j * p + b] * V_mat[j * p + b];
    }
    IC_bar[b] = IC_bar[b] / (1.0 * n_IC);
    IC_std[b] = sqrt(IC_std[b] / (1.0 * n_IC) - IC_bar[b] * IC_bar[b]);
  }

  free(V_mat);

}

void ALO_phase2(int *t_in, int *p_in, double *U_prev, double *X_now, double *sig2_in, double *IC_bar, double *IC_std, double *lam_in, double *U_now, double *W,
		int *q_in, double *Q) {
  double lam = lam_in[0];
  int p = p_in[0];
  int t = t_in[0];
  int q = q_in[0];
  int b;
  double *Ws;


  for (b = 0; b < p; b++) {
    U_now[b] = lam * X_now[b] + (1.0 - lam) * U_prev[b];
  }

  ALO_phase1(p_in, sig2_in, U_now, W);

  for (b = 0; b < p; b++) {
    W[b] = (2.0 - lam)/(lam * (1.0 - pow(1.0 - lam, 2*t))) * W[b];
  }

  Ws = (double *)malloc(p * sizeof(double));
  for (b = 0; b < p; b++) {
    Ws[b] = (W[b] - IC_bar[b])/IC_std[b];
  }

  Q[0] = kmax(Ws, q);

  free(Ws);

}

void SCAD_phase2(int *t_in, int *p_in, double *U_prev, double *X_now, double *sig2_in, double *a_in, double *IC_bar, double *IC_std, double *lam_in, double *U_now, double *V,
		int *q_in, double *Q) {
  double lam = lam_in[0];
  int p = p_in[0];
  int t = t_in[0];
  int q = q_in[0];
  int b;
  double *Vs;


  for (b = 0; b < p; b++) {
    U_now[b] = lam * X_now[b] + (1.0 - lam) * U_prev[b];
  }

  SCAD_phase1(p_in, sig2_in, a_in, U_now, V);

  for (b = 0; b < p; b++) {
    V[b] = (2.0 - lam)/(lam * (1.0 - pow(1.0 - lam, 2*t))) * V[b];
  }

  Vs = (double *)malloc(p * sizeof(double));
  for (b = 0; b < p; b++) {
    Vs[b] = (V[b] - IC_bar[b])/IC_std[b];
  }

  Q[0] = kmax(Vs, q);

  free(Vs);

}

void ALO_RL0(int *n_phase2_in, int *p_in, double *sig2_in, double *IC_bar, double *IC_std, double *lam_in, int *q_in, double *ctrlim_in, int *RL0) {
  int n_phase2 = n_phase2_in[0];
  int p = p_in[0];
  double ctrlim = ctrlim_in[0];
  int t, b;
  double *U_prev, *U_now, *X, *W, *Q;


  RL0[0] = 0;
  t = 0;
  Q = (double *)malloc(1 * sizeof(double));
  Q[0] = 0.0;
  
  U_prev = (double *)malloc(p * sizeof(double));
  U_now = (double *)malloc(p * sizeof(double));
  X = (double *)malloc(p * sizeof(double));
  W = (double *)malloc(p * sizeof(double));

  /* initialize U0 */

  for (b = 0; b < p; b++) {
    U_prev[b] = 0.0;
  }

  /* generate phase II observations and compute charting statistics */
  
  GetRNGstate();

  while ((t < n_phase2) && (Q[0] <= ctrlim)) {

    t = t + 1;
    for (b = 0; b < p; b++) {
      X[b] = sqrt(sig2_in[0]) * rnorm(0.0, 1.0);
    }
    ALO_phase2(&t, p_in, U_prev, X, sig2_in, IC_bar, IC_std, lam_in, U_now, W, q_in, Q);

    if (t <= 5) {
      Rprintf("Q = %g \n", Q[0]);
    }

    for (b = 0; b < p; b++) {
      U_prev[b] = U_now[b];
    }
  }

  PutRNGstate();

  RL0[0] = t;

  free(U_prev);
  free(U_now);
  free(X);
  free(W);
  free(Q);
    
}

void ALO_Qs(int *n_phase2_in, int *p_in, double *sig2_in, double *IC_bar, double *IC_std, double *lam_in, int *q_in, double *Qs) {
  int n_phase2 = n_phase2_in[0];
  int p = p_in[0];
  int t, b, t1;
  double *U_prev, *U_now, *X, *W, *Q;

  t = 0;
  Q = (double *)malloc(1 * sizeof(double));
  Q[0] = 0.0;
  
  U_prev = (double *)malloc(p * sizeof(double));
  U_now = (double *)malloc(p * sizeof(double));
  X = (double *)malloc(p * sizeof(double));
  W = (double *)malloc(p * sizeof(double));

  /* initialize U0 */

  for (b = 0; b < p; b++) {
    U_prev[b] = 0.0;
  }

  /* generate phase II observations and compute charting statistics */
  
  GetRNGstate();

  for (t = 0; t < n_phase2; t++) {

    t1 = t + 1; 
    for (b = 0; b < p; b++) {
      X[b] = sqrt(sig2_in[0]) * rnorm(0.0, 1.0);
    }
    ALO_phase2(&t1, p_in, U_prev, X, sig2_in, IC_bar, IC_std, lam_in, U_now, W, q_in, Q);
    Qs[t] = Q[0];
    /* if (t <= 5) { */
    /*   Rprintf("Q = %g \n", Q[0]); */
    /* } */

    for (b = 0; b < p; b++) {
      U_prev[b] = U_now[b];
    }
  }

  PutRNGstate();

  free(U_prev);
  free(U_now);
  free(X);
  free(W);
  free(Q);
    
}



void SCAD_RL0(int *n_phase2_in, int *p_in, double *sig2_in, double *a_in, double *IC_bar, double *IC_std, double *lam_in, int *q_in, double *ctrlim_in, int *RL0) {
  int n_phase2 = n_phase2_in[0];
  int p = p_in[0];
  double ctrlim = ctrlim_in[0];
  int t, b;
  double *U_prev, *U_now, *X, *V, *Q;


  RL0[0] = 0;
  t = 0;
  Q = (double *)malloc(1 * sizeof(double));
  Q[0] = 0.0;
  
  U_prev = (double *)malloc(p * sizeof(double));
  U_now = (double *)malloc(p * sizeof(double));
  X = (double *)malloc(p * sizeof(double));
  V = (double *)malloc(p * sizeof(double));

  /* initialize U0 */

  for (b = 0; b < p; b++) {
    U_prev[b] = 0.0;
  }

  /* generate phase II observations and compute charting statistics */
  
  GetRNGstate();

  while ((t < n_phase2) && (Q[0] <= ctrlim)) {

    t = t + 1;
    for (b = 0; b < p; b++) {
      X[b] = sqrt(sig2_in[0]) * rnorm(0.0, 1.0);
    }
    SCAD_phase2(&t, p_in, U_prev, X, sig2_in, a_in, IC_bar, IC_std, lam_in, U_now, V, q_in, Q);

    if (t <= 5) {
      Rprintf("Q = %g \n", Q[0]);
    }

    for (b = 0; b < p; b++) {
      U_prev[b] = U_now[b];
    }
  }

  PutRNGstate();

  RL0[0] = t;

  free(U_prev);
  free(U_now);
  free(X);
  free(V);
  free(Q);
    
}

void SCAD_Qs(int *n_phase2_in, int *p_in, double *sig2_in, double *a_in, double *IC_bar, double *IC_std, double *lam_in, int *q_in, double *Qs) {
  int n_phase2 = n_phase2_in[0];
  int p = p_in[0];
  int t, b, t1;
  double *U_prev, *U_now, *X, *V, *Q;

  t = 0;
  Q = (double *)malloc(1 * sizeof(double));
  Q[0] = 0.0;
  
  U_prev = (double *)malloc(p * sizeof(double));
  U_now = (double *)malloc(p * sizeof(double));
  X = (double *)malloc(p * sizeof(double));
  V = (double *)malloc(p * sizeof(double));

  /* initialize U0 */

  for (b = 0; b < p; b++) {
    U_prev[b] = 0.0;
  }

  /* generate phase II observations and compute charting statistics */
  
  GetRNGstate();

  for (t = 0; t < n_phase2; t++) {

    t1 = t + 1;
    for (b = 0; b < p; b++) {
      X[b] = sqrt(sig2_in[0]) * rnorm(0.0, 1.0);
    }
    SCAD_phase2(&t1, p_in, U_prev, X, sig2_in, a_in, IC_bar, IC_std, lam_in, U_now, V, q_in, Q);
    Qs[t] = Q[0];
    /* if (t <= 5) { */
    /*   Rprintf("Q = %g \n", Q[0]); */
    /* } */

    for (b = 0; b < p; b++) {
      U_prev[b] = U_now[b];
    }
  }

  PutRNGstate();

  free(U_prev);
  free(U_now);
  free(X);
  free(V);
  free(Q);
    
}


void ALO_ARL0(int *nrep_in, int *n_phase2_in, int *p_in, double *sig2_in, double *IC_bar, double *IC_std, double *lam_in, int *q_in, double *ctrlim_in, double *ARL0) {
  int nrep = nrep_in[0];
  int irep;
  int *RL0;

  RL0 = (int *)malloc(1 * sizeof(int));

  ARL0[0] = 0.0;

  GetRNGstate();

  for (irep = 0; irep < nrep; irep++) {

    ALO_RL0(n_phase2_in, p_in, sig2_in, IC_bar, IC_std, lam_in, q_in, ctrlim_in, RL0);

    /* Rprintf("RL[%i] = %i \n", irep, RL0[0]); */

    ARL0[0] = ARL0[0] + RL0[0];

  }

  PutRNGstate();

  ARL0[0] = ARL0[0] / (1.0 * nrep);

  free(RL0);
}

void ALO_Qmat(int *nrep_in, int *n_phase2_in, int *p_in, double *sig2_in, double *IC_bar, double *IC_std, double *lam_in, int *q_in, double *Qmat) {
  int nrep = nrep_in[0];
  int n_phase2 = n_phase2_in[0];
  int irep, t;
  double *Qs;

  Qs = (double *)malloc(n_phase2 * sizeof(double));

  GetRNGstate();

  for (irep = 0; irep < nrep; irep++) {

    ALO_Qs(n_phase2_in, p_in, sig2_in, IC_bar, IC_std, lam_in, q_in, Qs);
    for (t = 0; t < n_phase2; t++) {
      Qmat[irep * n_phase2 + t] = Qs[t];
    }

  }

  PutRNGstate();

  free(Qs);
}

void SCAD_ARL0(int *nrep_in, int *n_phase2_in, int *p_in, double *sig2_in, double *a_in, double *IC_bar, double *IC_std, double *lam_in, int *q_in, double *ctrlim_in,
	       double *ARL0) {
  int nrep = nrep_in[0];
  int irep;
  int *RL0;

  RL0 = (int *)malloc(1 * sizeof(int));

  ARL0[0] = 0.0;

  GetRNGstate();

  for (irep = 0; irep < nrep; irep++) {

    SCAD_RL0(n_phase2_in, p_in, sig2_in, a_in, IC_bar, IC_std, lam_in, q_in, ctrlim_in, RL0);

    /* Rprintf("RL[%i] = %i \n", irep, RL0[0]); */

    ARL0[0] = ARL0[0] + RL0[0];

  }

  PutRNGstate();

  ARL0[0] = ARL0[0] / (1.0 * nrep);

  free(RL0);
}

void SCAD_Qmat(int *nrep_in, int *n_phase2_in, int *p_in, double *sig2_in, double *a_in, double *IC_bar, double *IC_std, double *lam_in, int *q_in, double *Qmat) {
  int nrep = nrep_in[0];
  int n_phase2 = n_phase2_in[0];
  int irep, t;
  double *Qs;

  Qs = (double *)malloc(n_phase2 * sizeof(double));


  GetRNGstate();

  for (irep = 0; irep < nrep; irep++) {

    SCAD_Qs(n_phase2_in, p_in, sig2_in, a_in, IC_bar, IC_std, lam_in, q_in, Qs);
    for (t = 0; t < n_phase2; t++) {
      Qmat[irep * n_phase2 + t] = Qs[t];
    }

  }

  PutRNGstate();

  free(Qs);
}

void ALO_im(int *t_in, int *n1_in, int *n2_in, double *U_prev, double *Z_now, double *wcoef_mean, double *sig2_in, double *IC_bar, double *IC_std, double *lam_in,
	    int *q_in, double *U_now, double *Q) {
  int *J1_in, *J2_in;
  int NC1, NC2, b, p;
  double *W, *wcoef_now;

  J1_in = (int *)malloc(1 * sizeof(int));
  J2_in = (int *)malloc(1 * sizeof(int));

  J1_in[0] = (int)floor(log(n1_in[0] * 1.0) / log(2.0));
  J2_in[0] = (int)floor(log(n2_in[0] * 1.0) / log(2.0));

  NC1 = (int)round(pow(2, J1_in[0]));
  NC2 = (int)round(pow(2, J2_in[0]));

  p = NC1 * NC2;

  /* compute the wavelet coefficients for the current image */

  wcoef_now = (double *)malloc(p * sizeof(double));
  
  wav2d_coef(n1_in, n2_in, Z_now, J1_in, J2_in, wcoef_now);

  /* subtract true the wavelet coefficient mean */

  for (b = 0; b < p; b++) {
    wcoef_now[b] = wcoef_now[b] - wcoef_mean[b];
  }

  /* compute the LEWMA charting statistic on the centered wavelet coefficients */

  W = (double *)malloc(p * sizeof(double));

  ALO_phase2(t_in, &p, U_prev, wcoef_now, sig2_in, IC_bar, IC_std, lam_in, U_now, W, q_in, Q);

  free(J1_in);
  free(J2_in);
  free(wcoef_now);
  free(W);

}

void SCAD_im(int *t_in, int *n1_in, int *n2_in, double *U_prev, double *Z_now, double *wcoef_mean, double *sig2_in, double *a_in, double *IC_bar, double *IC_std, double *lam_in,
	    int *q_in, double *U_now, double *Q) {
  int *J1_in, *J2_in;
  int NC1, NC2, b, p;
  double *V, *wcoef_now;

  J1_in = (int *)malloc(1 * sizeof(int));
  J2_in = (int *)malloc(1 * sizeof(int));

  J1_in[0] = (int)floor(log(n1_in[0] * 1.0) / log(2.0));
  J2_in[0] = (int)floor(log(n2_in[0] * 1.0) / log(2.0));

  NC1 = (int)round(pow(2, J1_in[0]));
  NC2 = (int)round(pow(2, J2_in[0]));

  p = NC1 * NC2;

  /* compute the wavelet coefficients for the current image */

  wcoef_now = (double *)malloc(p * sizeof(double));
  
  wav2d_coef(n1_in, n2_in, Z_now, J1_in, J2_in, wcoef_now);

  /* subtract true the wavelet coefficient mean */

  for (b = 0; b < p; b++) {
    wcoef_now[b] = wcoef_now[b] - wcoef_mean[b];
  }

  /* compute the SCAD charting statistic on the centered wavelet coefficients */

  V = (double *)malloc(p * sizeof(double));

  SCAD_phase2(t_in, &p, U_prev, wcoef_now, sig2_in, a_in, IC_bar, IC_std, lam_in, U_now, V, q_in, Q);

  free(J1_in);
  free(J2_in);
  free(wcoef_now);
  free(V);

}

void change_loc(int *n_im_in, int *n1_in, int *n2_in, double *Z_seq, double *wcoef_mean, double *like, int *loc_hat) {
  int n_im = n_im_in[0];
  int n1 = n1_in[0];
  int n2 = n2_in[0];
  int j, b, p, i1, i2;
  int *J1_in, *J2_in;
  double *Z, *wcoef, *ybwd_sum;

  J1_in = (int *)malloc(1 * sizeof(int));
  J2_in = (int *)malloc(1 * sizeof(int));

  J1_in[0] = (int)floor(log(n1_in[0] * 1.0) / log(2.0));
  J2_in[0] = (int)floor(log(n2_in[0] * 1.0) / log(2.0));

  p = (int)round(pow(2, J1_in[0] + J2_in[0]));

  /* Rprintf("J1 = %i, J2 = %i, p = %i \n", J1_in[0], J2_in[0], p); */

  Z = (double *)malloc(n1 * n2 * sizeof(double));
  wcoef = (double *)malloc(p * sizeof(double));
  ybwd_sum = (double *)malloc(p * sizeof(double));

  /* initialize */
  
  for (b = 0; b < p; b++) {
    ybwd_sum[b] = 0.0;
  }

  /* backward looking */
  
  for (j = (n_im - 1); j >= 0; j--) {

    for (i1 = 0; i1 < n1; i1++) {
      for (i2 = 0; i2 < n2; i2++) {

	Z[i1 * n1 + i2] = Z_seq[j * n1 * n2 + i1 * n2 + i2];

      }
    }

    wav2d_coef(n1_in, n2_in, Z, J1_in, J2_in, wcoef);

    for (b = 0; b < p; b++) {
      ybwd_sum[b] = ybwd_sum[b] + wcoef[b] - wcoef_mean[b];
    }

    like[j] = 0.0;

    for (b = 0; b < p; b++) {
      like[j] = like[j] + ybwd_sum[b] * ybwd_sum[b];
    }

    like[j] = like[j] / (1.0 * (n_im - j)); 

  }

  /* change location MLE */

  loc_hat[0] = maxloc(n_im, like);

  free(J1_in);
  free(J2_in);
  free(Z);
  free(wcoef);
  free(ybwd_sum);

}

void ALO_solpath(int *p_in, double *y, int *q_in, double *mu_hat) {
  int p = p_in[0];
  int q = q_in[0];
  double *yabs;
  int b;
  double thresh;
  double const small_number = pow(10, -37);

  /* sort the absolute values */
  
  yabs = (double *)malloc(p * sizeof(double));
  for (b = 0; b < p; b++) {
    yabs[b] = fabs(y[b]);
  }

  quicksort(yabs, 0, p - 1);

  if (yabs[0] < small_number) {
    error("y contains a 0 value. \n");
  }

  for (b = 0; b < (p - 1); b++) {
    if ((yabs[b + 1] - yabs[b]) < small_number) {
      error("|y| contains ties. \n");
    }
  }

  if (q >= p) {
    thresh = 0.0;
  } else {
    thresh = pow(yabs[p - q - 1], 2);
  }

  for (b = 0; b < p; b++) {
    mu_hat[b] = thresholding(y[b] * y[b], thresh) / fabs(y[b]);
  }

  free(yabs);
}

void SCAD_solpath(int *p_in, double *y, int *q_in, double *a_in, double *mu_hat) {
  int p = p_in[0];
  int q = q_in[0];
  double a = a_in[0];
  double *yabs;
  int b;
  double thresh;
  double const small_number = pow(10, -37);

  /* sort the absolute values */
  
  yabs = (double *)malloc(p * sizeof(double));
  for (b = 0; b < p; b++) {
    yabs[b] = fabs(y[b]);
  }

  quicksort(yabs, 0, p - 1);

  if (yabs[0] < small_number) {
    error("y contains a 0 value. \n");
  }

  for (b = 0; b < (p - 1); b++) {
    if ((yabs[b + 1] - yabs[b]) < small_number) {
      error("|y| contains ties. \n");
    }
  }

  if (q >= p) {
    thresh = 0.0;
  } else {
    thresh = yabs[p - q - 1];
  }

  for (b = 0; b < p; b++) {
    mu_hat[b] = SCAD_1dmu(1.0, y[b], thresh, a);
  }

  free(yabs);
}

void ALO_change_wcoef(int *n_OCim_in, int *n1_in, int *n2_in, double *sig2_in, double *Z_seq, double *wcoef_mean, double *EBIC, double *OC_wcoef) {
  int n_OCim = n_OCim_in[0];
  int n1 = n1_in[0];
  int n2 = n2_in[0];
  double sig2 = sig2_in[0];
  int j, b, b1, p, i1, i2, q, q_min;
  int *J1_in, *J2_in;
  double *Z, *wcoef, *y_OCbar, *ALO_sol;
  double cc;

  J1_in = (int *)malloc(1 * sizeof(int));
  J2_in = (int *)malloc(1 * sizeof(int));

  J1_in[0] = (int)floor(log(n1_in[0] * 1.0) / log(2.0));
  J2_in[0] = (int)floor(log(n2_in[0] * 1.0) / log(2.0));

  p = (int)round(pow(2, J1_in[0] + J2_in[0]));

  Z = (double *)malloc(n1 * n2 * sizeof(double));
  wcoef = (double *)malloc(p * sizeof(double));
  y_OCbar = (double *)malloc(p * sizeof(double));

  /* initialize */
  
  for (b = 0; b < p; b++) {
    y_OCbar[b] = 0.0;
  }

  /* compute OC y average */
  
  for (j = 0; j < n_OCim; j++) {
    
    for (i1 = 0; i1 < n1; i1++) {
      for (i2 = 0; i2 < n2; i2++) {

	Z[i1 * n1 + i2] = Z_seq[j * n1 * n2 + i1 * n1 + i2];

      }
    }

    wav2d_coef(n1_in, n2_in, Z, J1_in, J2_in, wcoef);

    for (b = 0; b < p; b++) {
      y_OCbar[b] = y_OCbar[b] + wcoef[b] - wcoef_mean[b];
    }

  }

  for (b = 0; b < p; b++) {
    y_OCbar[b] = y_OCbar[b] / (1.0 * n_OCim);
  }

  /* compute adaptive lasso solution and the corresponding EBIC */

  ALO_sol = (double *)malloc(p * sizeof(double));

  cc = n_OCim / sig2;

  for (b = 0; b < p; b++) {
    q = b + 1;
    ALO_solpath(&p, y_OCbar, &q, ALO_sol);

    EBIC[b] = (2 * log(1.0 * p) + log(1.0 * n_OCim)) * q;
    for (b1 = 0; b1 < p; b1++) {
      EBIC[b] = EBIC[b] + cc * pow(y_OCbar[b1] - ALO_sol[b1], 2);
    }
    
  }

  /* find the OC wavelet coefficients by minimizing EBIC */

  q_min = minloc(p, EBIC) + 1;
  ALO_solpath(&p, y_OCbar, &q_min, OC_wcoef);

  free(J1_in);
  free(J2_in);
  free(Z);
  free(wcoef);
  free(y_OCbar);
  free(ALO_sol);
}

void SCAD_change_wcoef(int *n_OCim_in, int *n1_in, int *n2_in, double *sig2_in, double *a_in, double *Z_seq, double *wcoef_mean, double *EBIC,
		       double *OC_wcoef) {
  int n_OCim = n_OCim_in[0];
  int n1 = n1_in[0];
  int n2 = n2_in[0];
  double sig2 = sig2_in[0];
  int j, b, b1, p, i1, i2, q, q_min;
  int *J1_in, *J2_in;
  double *Z, *wcoef, *y_OCbar, *SCAD_sol;
  double cc;

  J1_in = (int *)malloc(1 * sizeof(int));
  J2_in = (int *)malloc(1 * sizeof(int));

  J1_in[0] = (int)floor(log(n1_in[0] * 1.0) / log(2.0));
  J2_in[0] = (int)floor(log(n2_in[0] * 1.0) / log(2.0));

  p = (int)round(pow(2, J1_in[0] + J2_in[0]));

  Z = (double *)malloc(n1 * n2 * sizeof(double));
  wcoef = (double *)malloc(p * sizeof(double));
  y_OCbar = (double *)malloc(p * sizeof(double));

  /* initialize */
  
  for (b = 0; b < p; b++) {
    y_OCbar[b] = 0.0;
  }

  /* compute OC y average */
  
  for (j = 0; j < n_OCim; j++) {
    
    for (i1 = 0; i1 < n1; i1++) {
      for (i2 = 0; i2 < n2; i2++) {

	Z[i1 * n1 + i2] = Z_seq[j * n1 * n2 + i1 * n1 + i2];

      }
    }

    wav2d_coef(n1_in, n2_in, Z, J1_in, J2_in, wcoef);

    for (b = 0; b < p; b++) {
      y_OCbar[b] = y_OCbar[b] + wcoef[b] - wcoef_mean[b];
    }

  }

  for (b = 0; b < p; b++) {
    y_OCbar[b] = y_OCbar[b] / (1.0 * n_OCim);
  }

  /* compute SCAD solution and the corresponding EBIC */

  SCAD_sol = (double *)malloc(p * sizeof(double));

  cc = n_OCim / sig2;

  for (b = 0; b < p; b++) {
    q = b + 1;
    SCAD_solpath(&p, y_OCbar, &q, a_in, SCAD_sol);

    EBIC[b] = (2 * log(1.0 * p) + log(1.0 * n_OCim)) * q;
    for (b1 = 0; b1 < p; b1++) {
      EBIC[b] = EBIC[b] + cc * pow(y_OCbar[b1] - SCAD_sol[b1], 2);
    }
    
  }

  /* find the OC wavelet coefficients by minimizing EBIC */

  q_min = minloc(p, EBIC) + 1;
  SCAD_solpath(&p, y_OCbar, &q_min, a_in, OC_wcoef);

  free(J1_in);
  free(J2_in);
  free(Z);
  free(wcoef);
  free(y_OCbar);
  free(SCAD_sol);
}

void w2d_indx(int *J1_in, int *J2_in, int *b_in, int *p_2d, int *j_2d, int *k_2d) {
  int J1 = J1_in[0];
  int J2 = J2_in[0];
  int b = b_in[0];
  int NC1, NC2, di, pp;

  NC1 = (int)round(pow(2, J1));
  NC2 = (int)round(pow(2, J2));

  if (b > (NC1 * NC2)) {
    error("b > pow(2, J1 + J2). \n");
  }

  p_2d[0] = b / NC2;
  p_2d[1] = b % NC2;

  for (di = 0; di < 2; di++) {
    pp = p_2d[di];
    if (pp == 0) {
      j_2d[di] = 0;
      k_2d[di] = 0;
    } else {
      j_2d[di] = (int)(floor(log(1.0 * pp)/log(2.0)));
      k_2d[di] = pp - (int)(round(pow(2, j_2d[di])));
    }
  }

}

void w2d_pxl(int *n1_in, int *n2_in, int *J1_in, int *J2_in, int *nb_in, int *bseq, int *pxl) {
  int n1 = n1_in[0];
  int n2 = n2_in[0];
  int nb = nb_in[0];
  int ib, b, i1, i2;
  int *p_2d, *j_2d, *k_2d;
  double x, y;

  /* initialize the pixel matrix to be all 0 */

  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {

      pxl[i1 * n2 + i2] = 0;

    }
  }

  /* convert wavelet coefficients to pixel locations */
  
  p_2d = (int *)malloc(2 * sizeof(int));
  j_2d = (int *)malloc(2 * sizeof(int));
  k_2d = (int *)malloc(2 * sizeof(int));

  for (ib = 0; ib < nb; ib++) {
    b = bseq[ib];
    w2d_indx(J1_in, J2_in, &b, p_2d, j_2d, k_2d);

    for (i1 = 0; i1 < n1; i1++) {
      x = (1.0 * i1) / n1;
      for (i2 = 0; i2 < n2; i2++) {
	y = (1.0 * i2) / n2;

	if ((x >= k_2d[0]/pow(2.0, j_2d[0])) && (x < (k_2d[0] + 1)/pow(2.0, j_2d[0])) && (y >= k_2d[1]/pow(2.0, j_2d[1])) && (y < (k_2d[1] + 1)/pow(2.0, j_2d[1]))) {
	  pxl[i1 * n2 + i2] = 1;
	}

      }
    }

  }

  free(p_2d);
  free(j_2d);
  free(k_2d);

}

void dQ(int *n1_in, int *n2_in, int *D, int *Dhat, double *w_in, double *dq) {
  int n1 = n1_in[0];
  int n2 = n2_in[0];
  double w = w_in[0];
  int i1, i2, TN, TP, FP, FN;

  /* true positives */
  
  TP = 0;
  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {
      TP = TP + D[i1 * n2 + i2];
    }
  }

  /* true negatives */

  TN = n1 * n2  - TP;

  /* false positives */

  FP = 0;
  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {

      if ((Dhat[i1 * n2 + i2] == 1) && (D[i1 * n2 + i2] != 1)) {
	FP = FP + 1;
      }

    }
  }

  /* false negatives */

  FN = 0;
  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {
      
      if ((D[i1 * n2 + i2] == 1) && (Dhat[i1 * n2 + i2] != 1)) {
	FN = FN + 1;
      }

    }
  }

  /* compute dq = w FP/TN + (1-w) FN/TP */

  if ((TP >= 1) && (TN >= 1)) {
    dq[0] = w * (1.0 * FP) / TN + (1.0 - w) * (1.0 * FN) / TP;
  } else {
    if (TP == 0) {
      dq[0] = (1.0 * FP) / TN;
    } else { 			/* TN = 0 */
      dq[0] = (1.0 * FN) / TP;
    }
  }

}

void w1d_indx(int *J_in, int *b_in, int *j, int *k) {
  int J = J_in[0];
  int b = b_in[0];
  int NC;

  NC = (int)round(pow(2, J));
  if (b > NC) {
    error("b > pow(2, J). \n");
  }

  if (b == 0) {
    j[0] = 0;
    k[0] = 0;
  } else {
    j[0] = (int)(floor(log(1.0 * b)/log(2.0)));
    k[0] = b - (int)(round(pow(2, j[0])));
  }
  

}

void w1ds_pxl(int *n1_in, int *n2_in, int *J2_in, int *BETA01, int *pxl) {
  int n1 = n1_in[0];
  int n2 = n2_in[0];
  int J2 = J2_in[0];
  int nb, ib, b, i1, i2;
  int *j, *k;
  double y;

  /* initialize the pixel vector to be all 0 */

  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {
      pxl[i1 * n2 + i2] = 0;
    }
  }

  /* convert 1-d wavelet coefficients to pixel locations */

  nb = (int)(pow(2, J2));
  
  j = (int *)malloc(1 * sizeof(int));
  k = (int *)malloc(1 * sizeof(int));

  for (i1 = 0; i1 < n1; i1++) {
    for (ib = 0; ib < nb; ib++) {

      if (BETA01[i1 * nb + ib] == 1) {
	b = ib;
	w1d_indx(J2_in, &b, j, k);

	for (i2 = 0; i2 < n2; i2++) {
	  y = (1.0 * i2) / n2;

	  if ((y >= k[0]/pow(2.0, j[0])) && (y < (k[0] + 1)/pow(2.0, j[0]))) {
	    pxl[i1 * n2 + i2] = 1;
	  }

	}

      }

    }
  }

  free(j);
  free(k);

}

double max(int n, double *x) {
  int i;
  double xmax;

  xmax = x[0];
  for (i = 0; i < n; i++) {
    if (xmax < x[i]) {
      xmax = x[i];
    }
  }
  return(xmax);
}

void maxarg(int n, double *x, double *maxval, int *maxloc) {
  int i;

  maxval[0] = x[0];
  maxloc[0] = 0;
  for (i = 0; i < n; i++) {
    if (maxval[0] < x[i]) {
      maxval[0] = x[i];
      maxloc[0] = i;
    }
  }
}


void GLR1d_RLs0(int *nRuns_in, int *upbnd_in, int *m_in, double *UCL_in, double *RLs0) {
  int nRuns = nRuns_in[0];
  int upbnd = upbnd_in[0];
  int m = m_in[0];
  double UCL = UCL_in[0];
  int s, t, iRun;
  double *Rms, *bwdsum;
  double z, maxRms;

  Rms = (double *)malloc(m * sizeof(double));
  bwdsum = (double *)malloc(upbnd * sizeof(double));

  GetRNGstate();

  for (iRun = 0; iRun < nRuns; iRun++) {

    for (t = 0; t < m; t++) {
      Rms[t] = 0.0;
    }
    for (t = 0; t < upbnd; t++) {
      bwdsum[t] = 0.0;
    }

    s = 0;
    maxRms = 0.0;

    while (s < upbnd && maxRms < UCL) {

      z = rnorm(0.0, 1.0);
      /* Rprintf("z[%i] = %g \n", s, z); */

      for (t = 0; t <= s; t++) {
	bwdsum[t] = bwdsum[t] + z;
      }

      if (s < m) {
	for (t = 0; t <= s; t++) {
	  Rms[t] = 0.5 * (s + 1 - t) * pow(bwdsum[t] / (1.0 * (s + 1 - t)), 2);
	}
      } else {
	for (t = 0; t < m; t++) {
	  Rms[t] = 0.5 * (m - t) * pow(bwdsum[s - (m - 1 - t)] / (1.0 * (m - t)), 2);
	}
      }

      maxRms = max(m, Rms);
      /* Rprintf("max Rms = %g \n", maxRms); */
    
      s = s + 1;
    }

    RLs0[iRun] = 1.0 * s;
  }

  PutRNGstate();

  free(Rms);
  free(bwdsum);

}


void IdpGLR_RLs0(int *nGLRs_in, int *nRuns_in, int *upbnd_in, int *m_in, double *UCL_in, double *RLs0) {
  int nGLRs = nGLRs_in[0];
  int nRuns = nRuns_in[0];
  int upbnd = upbnd_in[0];
  int m = m_in[0];
  double UCL = UCL_in[0];
  int s, t, iRun, iGLR;
  double *Rms, *bwdsum, *maxRms, *tempRms;
  double z, maxmaxRms;

  Rms = (double *)malloc(m * nGLRs * sizeof(double));
  bwdsum = (double *)malloc(upbnd * nGLRs * sizeof(double));
  maxRms = (double *)malloc(nGLRs * sizeof(double));
  tempRms = (double *)malloc(m * sizeof(double));

  GetRNGstate();

  for (iRun = 0; iRun < nRuns; iRun++) {

    for (iGLR = 0; iGLR < nGLRs; iGLR++) {
      maxRms[iGLR] = 0.0;
      for (t = 0; t < m; t++) {
	Rms[iGLR * m + t] = 0.0;
      }
      for (t = 0; t < upbnd; t++) {
	bwdsum[iGLR * upbnd + t] = 0.0;
      } 
    }

    s = 0;
    maxmaxRms = 0;

    while (s < upbnd && maxmaxRms < UCL) {

      for (iGLR = 0; iGLR < nGLRs; iGLR++) {
	
	z = rnorm(0.0, 1.0);
	/* Rprintf("z[%i] = %g \n", s, z); */

	for (t = 0; t <= s; t++) {
	  bwdsum[iGLR * upbnd + t] = bwdsum[iGLR * upbnd + t] + z;
	}

	if (s < m) {
	  for (t = 0; t <= s; t++) {
	    Rms[iGLR * m + t] = 0.5 * (s + 1 - t) * pow(bwdsum[iGLR * upbnd + t] / (1.0 * (s + 1 - t)), 2);
	  }
	} else {
	  for (t = 0; t < m; t++) {
	    Rms[iGLR * m + t] = 0.5 * (m - t) * pow(bwdsum[iGLR * upbnd + s - (m - 1 - t)] / (1.0 * (m - t)), 2);
	  }
	}

	for (t = 0; t < m; t++) {
	  tempRms[t] = Rms[iGLR * m + t];
	}

	maxRms[iGLR] = max(m, tempRms);
	/* Rprintf("max Rms = %g \n", maxRms[iGLR]); */
      }

      maxmaxRms = max(nGLRs, maxRms);
      /* Rprintf("max maxRms = %g \n", maxmaxRms); */
      s = s + 1;
    }

    RLs0[iRun] = 1.0 * s;
  }

  PutRNGstate();

  free(Rms);
  free(bwdsum);
  free(maxRms);
  free(tempRms);

}

void GLRs_phase2(int *nGLRs_in, int *m_in, int *s_in, double *sig2_in, double *BSprev, double *Xnow, double *Xmean, double *BSnow, double *GLRstats, int *GLRtaus) {
  int nGLRs = nGLRs_in[0];
  int m = m_in[0];
  int s = s_in[0];
  double sig2 = sig2_in[0];
  int iGLR, t;
  double *RR, *RRmaxval;
  int *RRmaxloc;

  /* update the backward sum matrix */

  if (s <= m) {
    for (t = 0; t < s; t++) {
      for (iGLR = 0; iGLR < nGLRs; iGLR++) {
	BSnow[iGLR * m + t] = BSprev[iGLR * m + t] + Xnow[iGLR];
      }
    }
  } else {			/* s >= m  */
    for (t = 0; t < m; t++) {
      for (iGLR = 0; iGLR < nGLRs; iGLR++) {
	if (t < (m - 1)) {
	  BSnow[iGLR * m + t] = BSprev[iGLR * m + (t + 1)] + Xnow[iGLR];
	} else {
	  BSnow[iGLR * m + t] = Xnow[iGLR];
	}
      }
    }
  }

  /* compute all the GLR statistics */

  RR = (double *)malloc(m * sizeof(double));
  for (t = 0; t < m; t++) {
    RR[t] = 0.0;
  }
  RRmaxval = (double *)malloc(1 * sizeof(double));
  RRmaxloc = (int *)malloc(1 * sizeof(int));
  
  for (iGLR = 0; iGLR < nGLRs; iGLR++) {
    if (s <= m) {
      for (t = 0; t < s; t++) {
	RR[t] = (0.5 * (s - t) / sig2) * pow(BSnow[iGLR * m + t] / (1.0 * (s - t)) - Xmean[iGLR], 2);
      }
    } else {
      for (t = 0; t < m; t++) {
	RR[t] = (0.5 * (m - t) / sig2) * pow(BSnow[iGLR * m + t] / (1.0 * (m - t)) - Xmean[iGLR], 2);
      }
    }

    maxarg(m, RR, RRmaxval, RRmaxloc);
    GLRstats[iGLR] = RRmaxval[0];
    GLRtaus[iGLR] = RRmaxloc[0];
  }

  free(RR);
  free(RRmaxval);
  free(RRmaxloc);
}
