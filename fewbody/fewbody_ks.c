/* -*- linux-c -*- */
/* fewbody_ks.c

   Copyright (C) 2002-2004 John M. Fregeau
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include "fewbody.h"

/* the dot product of two K-S vectors */
inline double fb_ks_dot(double x[4], double y[4])
{
	return(x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3]);
}

/* the modulus of a K-S vector */
inline double fb_ks_mod(double x[4])
{
	return(sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]));
}

/* calculate Q given q */
void fb_calc_Q(double q[4], double Q[4])
{
	double u[4];
	
	u[0] = sqrt(0.5 * (fb_mod(q) + fabs(q[0])));
	u[1] = 0.5 * q[1] / u[0];
	u[2] = 0.5 * q[2] / u[0];
	u[3] = 0.0;

	if (q[0] > 0.0) {
		Q[0] = u[0];
		Q[1] = u[1];
		Q[2] = u[2];
		Q[3] = u[3];
	} else {
		Q[0] = u[1];
		Q[1] = u[0];
		Q[2] = u[3];
		Q[3] = u[2];
	}
}

/* calculate the K-S matrix */
void fb_calc_ksmat(double Q[4], double Qmat[4][4])
{
	Qmat[0][0] = Q[0];
	Qmat[0][1] = -Q[1];
	Qmat[0][2] = -Q[2];
	Qmat[0][3] = Q[3];
	
	Qmat[1][0] = Q[1];
	Qmat[1][1] = Q[0];
	Qmat[1][2] = -Q[3];
	Qmat[1][3] = -Q[2];
	
	Qmat[2][0] = Q[2];
	Qmat[2][1] = Q[3];
	Qmat[2][2] = Q[0];
	Qmat[2][3] = Q[1];
	
	Qmat[3][0] = Q[3];
	Qmat[3][1] = -Q[2];
	Qmat[3][2] = Q[1];
	Qmat[3][3] = -Q[0];
}

/* calculate the a matrix */
void fb_calc_amat(double **a, int nstar, int kstar)
{
	int i, j, k;
	
	/* first zero out the matrix */
	for (i=0; i<nstar; i++) {
		for (k=0; k<kstar; k++) {
			a[i][k] = 0.0;
		}
	}
	
	/* then set the non-zero components */
	k = -1;
	for (i=0; i<nstar-1; i++) {
		for (j=i+1; j<nstar; j++) {
			k++;
			a[i][k] = 1.0;
			a[j][k] = -1.0;
		}
	}
}

/* calculate the T matrix */
void fb_calc_Tmat(double **a, double *m, double **T, int nstar, int kstar)
{
	int u, v, e;
	
	for (u=0; u<kstar; u++) {
		for (v=0; v<kstar; v++) {
			T[u][v] = 0.0;
			for (e=0; e<nstar; e++) {
				T[u][v] += a[e][u] * a[e][v] / m[e];
			}
			T[u][v] *= 0.5;
		}
	}
}

/* the derivatives function for the GSL ODE integrator */
int fb_ks_func(double s, const double *y, double *f, void *params)
{
	int k, l, m, kstar;
	double *M, **amat, **Tmat, Einit;
	double **Q, **P, **p, **A, **TP, **TQ, **UQ, *d;
	double Qmat[4][4], Pmat[4][4], Astar[4], T, U, val, L, H, G, GT, GU;

	/* set parameters */
	kstar = (*(fb_ks_params_t *) params).kstar;
	M = (*(fb_ks_params_t *) params).M;
	amat = (*(fb_ks_params_t *) params).amat;
	Tmat = (*(fb_ks_params_t *) params).Tmat;
	Einit = (*(fb_ks_params_t *) params).Einit;

	/* allocate memory */
	Q = fb_malloc_matrix(kstar, 4);
	P = fb_malloc_matrix(kstar, 4);
	p = fb_malloc_matrix(kstar, 4);
	A = fb_malloc_matrix(kstar, 4);
	TP = fb_malloc_matrix(kstar, 4);
	TQ = fb_malloc_matrix(kstar, 4);
	UQ = fb_malloc_matrix(kstar, 4);
	d = fb_malloc_vector(kstar);

	/* set Q_k, P_k, and p_k */
	for (k=0; k<kstar; k++) {
		/* read Q_k and P_k directly from y[] */
		for (l=0; l<4; l++) {
			Q[k][l] = y[k*8+1+l];
			P[k][l] = y[k*8+4+1+l];
		}
		
		/* calculate the KS matrix from Q_k */
		fb_calc_ksmat(Q[k], Qmat);
		
		/* and then calculate p_k */
		for (l=0; l<4; l++) {
			p[k][l] = 0.0;
			for (m=0; m<4; m++) {
				p[k][l] += Qmat[l][m] * P[k][m];
			}
			p[k][l] /= 2.0 * fb_ks_dot(Q[k], Q[k]);
		}
	}
	
	/* set intermediate variables */
	T = 0.0;
	U = 0.0;
	for (k=0; k<kstar; k++) {
		/* first calculate the K-S matrix for this value of k */
		fb_calc_ksmat(Q[k], Qmat);
		
		/* calculate A_k */
		for (m=0; m<4; m++) {
			A[k][m] = 0.0;
			for (l=0; l<kstar; l++) {
				A[k][m] += Tmat[k][l] * p[l][m];
			}
		}
		
		/* then d_k */
		d[k] = fb_ks_dot(A[k], p[k]);

		/* increment the sums for T and U */
		T += d[k];
		U += M[k] / fb_ks_dot(Q[k], Q[k]);
		
		/* then calculate the first partial derivative of T with respect to P_k */
		for (m=0; m<4; m++) {
			TP[k][m] = 0.0;
			for (l=0; l<4; l++) {
				TP[k][m] += Qmat[l][m] * A[k][l];
			}
			TP[k][m] /= fb_ks_dot(Q[k], Q[k]);
		}
		
		/* A^*_k */
		Astar[0] = A[k][0];
		Astar[1] = A[k][1];
		Astar[2] = A[k][2];
		Astar[3] = -A[k][3];
		
		/* calculate the K-S matrix for P_k */
		fb_calc_ksmat(P[k], Pmat);
		
		/* then calculate the first partial derivatives of T and U */
		for (l=0; l<4; l++) {
			val = 0.0;
			for (m=0; m<4; m++) {
				val += Pmat[m][l] * Astar[m];
			}
			TQ[k][l] = (val - 4.0 * d[k] * Q[k][l]) / fb_ks_dot(Q[k], Q[k]);
			UQ[k][l] = -2.0 * M[k] * Q[k][l] / fb_sqr(fb_ks_dot(Q[k], Q[k]));
		}
	}
	
	/* set the Lagrangian, Hamiltonian, and the initial energy */
	L = T + U;
	H = T - U;
	G = (H - Einit) / L;
	GT = (1.0 - G) / L;
	GU = -(1.0 + G) / L;
	
	/* set derivatives */
	f[0] = 1.0 / L;
	for (k=0; k<kstar; k++) {
		for (l=0; l<4; l++) {
			f[k*8+1+l] = GT * TP[k][l];
			f[k*8+4+1+l] = -GT * TQ[k][l] - GU * UQ[k][l];
		}
	}

	/* free memory */
	fb_free_matrix(Q);
	fb_free_matrix(P);
	fb_free_matrix(p);
	fb_free_matrix(A);
	fb_free_matrix(TP);
	fb_free_matrix(TQ);
	fb_free_matrix(UQ);
	fb_free_vector(d);

	/* all done */
	return(GSL_SUCCESS);
}

/* function to calculate the Einit parameter for the integrator */
/* the code here is a verbatim copy from the first part of fb_ks_func */
double fb_ks_Einit(const double *y, fb_ks_params_t params)
{
	int k, l, m;
	double **Q, **P, **p, **A, *d;
	double Qmat[4][4], T, U;

	/* allocate memory */
	Q = fb_malloc_matrix(params.kstar, 4);
	P = fb_malloc_matrix(params.kstar, 4);
	p = fb_malloc_matrix(params.kstar, 4);
	A = fb_malloc_matrix(params.kstar, 4);
	d = fb_malloc_vector(params.kstar);

	/* set Q_k, P_k, and p_k */
	for (k=0; k<params.kstar; k++) {
		/* read Q_k and P_k directly from y[] */
		for (l=0; l<4; l++) {
			Q[k][l] = y[k*8+1+l];
			P[k][l] = y[k*8+4+1+l];
		}
		
		/* calculate the KS matrix from Q_k */
		fb_calc_ksmat(Q[k], Qmat);
		
		/* and then calculate p_k */
		for (l=0; l<4; l++) {
			p[k][l] = 0.0;
			for (m=0; m<4; m++) {
				p[k][l] += Qmat[l][m] * P[k][m];
			}
			p[k][l] /= 2.0 * fb_ks_dot(Q[k], Q[k]);
		}
	}
	
	/* set intermediate variables */
	T = 0.0;
	U = 0.0;
	for (k=0; k<params.kstar; k++) {
		/* first calculate the K-S matrix for this value of k */
		fb_calc_ksmat(Q[k], Qmat);
		
		/* calculate A_k */
		for (m=0; m<4; m++) {
			A[k][m] = 0.0;
			for (l=0; l<params.kstar; l++) {
				A[k][m] += params.Tmat[k][l] * p[l][m];
			}
		}
		
		/* then d_k */
		d[k] = fb_ks_dot(A[k], p[k]);

		/* increment the sums for T and U */
		T += d[k];
		U += params.M[k] / fb_ks_dot(Q[k], Q[k]);
	}
	
	/* free memory */
	fb_free_matrix(Q);
	fb_free_matrix(P);
	fb_free_matrix(p);
	fb_free_matrix(A);
	fb_free_vector(d);
	
	return(T-U);
}

/* function to convert from Euclidean coordinates to K-S coordinates */
void fb_euclidean_to_ks(fb_obj_t **star, double *y, int nstar, int kstar)
{
	int i, j, k, l, m;
	double **q, **p, Q[4], P[4], Qmat[4][4];

	q = fb_malloc_matrix(kstar, 4);
	p = fb_malloc_matrix(kstar, 4);

	/* then calculate q_k and p_k */
	k = -1;
	for (i=0; i<nstar-1; i++) {
		for (j=i+1; j<nstar; j++) {
			k++;
			for (l=0; l<3; l++) {
				q[k][l] = star[i]->x[l] - star[j]->x[l];
				p[k][l] = (star[i]->m * star[i]->v[l] - star[j]->m * star[j]->v[l])/((double) nstar);
			}
			q[k][3] = 0.0;
			p[k][3] = 0.0;
		}
	}

	/* and then Q_k and P_k */
	for (k=0; k<kstar; k++) {
		/* calculate Q_k */
		fb_calc_Q(q[k], Q);

		/* then P_k */
		fb_calc_ksmat(Q, Qmat);
		for (l=0; l<4; l++) {
			P[l] = 0.0;
			for (m=0; m<4; m++) {
				P[l] += Qmat[m][l] * p[k][m];
			}
			P[l] *= 2.0;
			
			/* then set y_i */
			y[8*k+l+1] = Q[l];
			y[8*k+l+4+1] = P[l];
		}
	}

	fb_free_matrix(q);
	fb_free_matrix(p);
}

/* function to convert from K-S coordinates to Euclidean coordinates */
void fb_ks_to_euclidean(double *y, fb_obj_t **star, int nstar, int kstar)
{
	int i, j, k, l;
	double *m, mtot, **Q, **P, **q, **p, Qmat[4][4];

	/* allocate memory */
	m = fb_malloc_vector(nstar);
	Q = fb_malloc_matrix(kstar, 4);
	P = fb_malloc_matrix(kstar, 4);
	q = fb_malloc_matrix(kstar, 4);
	p = fb_malloc_matrix(kstar, 4);

	/* set the masses */
	mtot = 0.0;
	for (i=0; i<nstar; i++) {
		m[i] = star[i]->m;
		mtot += m[i];
	}

	/* set Q_k, P_k, and p_k */
	for (k=0; k<kstar; k++) {
		/* read Q_k and P_k directly from y[] */
		for (l=0; l<4; l++) {
			Q[k][l] = y[k*8+1+l];
			P[k][l] = y[k*8+4+1+l];
		}
		
		/* calculate the KS matrix from Q_k */
		fb_calc_ksmat(Q[k], Qmat);
		
		/* and then calculate q_k and p_k */
		for (l=0; l<4; l++) {
			q[k][l] = 0.0;
			p[k][l] = 0.0;
			for (j=0; j<4; j++) {
				q[k][l] += Qmat[l][j] * Q[k][j];
				p[k][l] += Qmat[l][j] * P[k][j];
			}
			p[k][l] /= 2.0 * fb_ks_dot(Q[k], Q[k]);
		}
	}
	
	/* set r_i and v_i */
	for (i=0; i<nstar; i++) {
		for (k=0; k<3; k++) {
			star[i]->x[k] = 0.0;
			star[i]->v[k] = 0.0;
		}
		
		for (j=i+1; j<nstar; j++) {
			for (k=0; k<3; k++) {
				star[i]->x[k] += m[j] * q[FB_KS_K(i, j, nstar)][k];
				star[i]->v[k] += p[FB_KS_K(i, j, nstar)][k];
			}
		}

		for (j=0; j<i; j++) {
			for (k=0; k<3; k++) {
				star[i]->x[k] += - m[j] * q[FB_KS_K(j, i, nstar)][k];
				star[i]->v[k] += - p[FB_KS_K(j, i, nstar)][k];
			}
		}

		for (k=0; k<3; k++) {
			star[i]->x[k] /= mtot;
			star[i]->v[k] /= m[i];
		}
	}
	
	/* free memory */
	fb_free_vector(m);
	fb_free_matrix(Q);
	fb_free_matrix(P);
	fb_free_matrix(q);
	fb_free_matrix(p);
}
