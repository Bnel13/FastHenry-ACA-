/*!\page LICENSE LICENSE

Copyright (C) 2003 by the Board of Trustees of Massachusetts Institute of
Technology, hereafter designated as the Copyright Owners.

License to use, copy, modify, sell and/or distribute this software and
its documentation for any purpose is hereby granted without royalty,
subject to the following terms and conditions:

1.  The above copyright notice and this permission notice must
appear in all copies of the software and related documentation.

2.  The names of the Copyright Owners may not be used in advertising or
publicity pertaining to distribution of the software without the specific,
prior written permission of the Copyright Owners.

3.  THE SOFTWARE IS PROVIDED "AS-IS" AND THE COPYRIGHT OWNERS MAKE NO
REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, BY WAY OF EXAMPLE, BUT NOT
LIMITATION.  THE COPYRIGHT OWNERS MAKE NO REPRESENTATIONS OR WARRANTIES OF
MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE
SOFTWARE WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS TRADEMARKS OR OTHER
RIGHTS. THE COPYRIGHT OWNERS SHALL NOT BE LIABLE FOR ANY LIABILITY OR DAMAGES
WITH RESPECT TO ANY CLAIM BY LICENSEE OR ANY THIRD PARTY ON ACCOUNT OF, OR
ARISING FROM THE LICENSE, OR ANY SUBLICENSE OR USE OF THE SOFTWARE OR ANY
SERVICE OR SUPPORT.

LICENSEE shall indemnify, hold harmless and defend the Copyright Owners and
their trustees, officers, employees, students and agents against any and all
claims arising out of the exercise of any rights under this Agreement,
including, without limiting the generality of the foregoing, against any
damages, losses or liabilities whatsoever with respect to death or injury to
person or damage to property arising from or out of the possession, use, or
operation of Software or Licensed Program(s) by LICENSEE or its customers.

*/



/* this file contains the functions for exact calculation of the
 self-inductance of a rectangular bar and the mutual inductance of
 two parallel rectangular bars.

 Also, it contains the code for the lookup table for these inductances */

#include "induct.h"

/* this counts on all arguments being positive! */
#define compare(x,y,eps) (  (((x)==0 && (y)==0) || (fabs((x) - (y)) < eps*((x) + (y)) )) \
  ? 0 : (  ((x) > (y)) ? 1 : -1 )  )

#define nearzero(x) (fabs(x) < EPS)

/* SRW */
double self(double, double, double);
int edges_parallel(FILAMENT*, FILAMENT*, double*, int*);
void get_wid(FILAMENT*, double*);
void get_height(FILAMENT*, double*, double*);
double exact_mutual(FILAMENT*, FILAMENT*, int, double*, double*,
		enum degen_type, enum degen_type);
void fill_4(double*, double, double, double);
double eval_eq(double, double, double, double);
double log_term(double, double, double, double, double);
double tan_term(double, double, double, double, double);
int lookup(FILAMENT*, FILAMENT*, int, double*, double*, double*, double*, int*,
		Table***, int*);
void fill_dims(FILAMENT*, FILAMENT*, double*, double*, double*, int);
void fill_dims_seg(FILAMENT*, FILAMENT*, double*, double*, double*, int);
int find_dims(double*, int, Table**, double*, int*, Table***);
void put_in_table(FILAMENT*, FILAMENT*, int, double, double*, int, Table**, int);
void init_table(void);
int get_table_mem(void);
void destroy_table(void);
char *AllocAnEntry(AllocInfo*);
void DestroyEntries(AllocInfo*);
int MemoryForEntries(AllocInfo*);
double brick_to_brick(double, double, double, double, double, double, double,
		double, double);
double flat_to_flat_tape(double, double, double, double, double, double, double);
double eval_eq_tape(double, double, double, double);
double flat_to_skinny_tape(double, double, double, double, double, double,
		double);
double eval_eq_tape2(double, double, double, double);
double tan_tape(double, double, double, double);
double tape_to_fil(double, double, double, double, double, double);
double brick_to_fil(double, double, double, double, double, double, double);

/* self inductance */
double self(double W, double L, double T) {

	double w, t, aw, at, ar, r, z;

	w = W / L;
	t = T / L;
	r = sqrt(w * w + t * t);
	aw = sqrt(w * w + 1.0);
	at = sqrt(t * t + 1.0);
	ar = sqrt(w * w + t * t + 1.0);

	z =
			0.25
					* ((1 / w) * asinh(w / at) + (1 / t) * asinh(t / aw)
							+ asinh(1 / r));
	z +=
			(1 / 24.0)
					* ((t * t / w) * asinh(w / (t * at * (r + ar)))
							+ (w * w / t) * asinh(t / (w * aw * (r + ar)))
							+ ((t * t) / (w * w))
									* asinh(w * w / (t * r * (at + ar)))
							+ ((w * w) / (t * t))
									* asinh(t * t / (w * r * (aw + ar)))
							+ (1.0 / (w * t * t))
									* asinh(w * t * t / (at * (aw + ar)))
							+ (1.0 / (t * w * w))
									* asinh(t * w * w / (aw * (at + ar))));
	z -= (1.0 / 6.0)
			* ((1.0 / (w * t)) * atan(w * t / ar) + (t / w) * atan(w / (t * ar))
					+ (w / t) * atan(t / (w * ar)));
	z -= (1.0 / 60.0)
			* (((ar + r + t + at) * t * t)
					/ ((ar + r) * (r + t) * (t + at) * (at + ar))
					+ ((ar + r + w + aw) * (w * w))
							/ ((ar + r) * (r + w) * (w + aw) * (aw + ar))
					+ (ar + aw + 1 + at)
							/ ((ar + aw) * (aw + 1) * (1 + at) * (at + ar)));
	z -= (1.0 / 20.0)
			* ((1.0 / (r + ar)) + (1.0 / (aw + ar)) + (1.0 / (at + ar)));

	z *= (2.0 / PI);
	z *= L; /* this is inductance */

	return z;

}

/* This assumes the lengths of the fils are parallel and returns 1 if the 
 side faces are parallel also
 */
int edges_parallel(FILAMENT *fil_j, FILAMENT *fil_m, double *wid1, int *whperp) {
	double *wid_j = fil_j->segm->widthdir;
	double *wid_m = fil_m->segm->widthdir;
	double wid2[3];
	double prod, mj, mm;

	if (wid_j == NULL && wid_m == NULL) {
		/* both have unspecified width directions and their lengths are assumed
		 parallel, so they are parallel */
		*whperp = 0;
		return 1 == 1;
	} else {
		if (wid_j == NULL) {
			/* get_wid(fil_j, wid1); */
			wid_j = wid1;
		}
		if (wid_m == NULL) {
			get_wid(fil_m, wid2);
			wid_m = wid2;
		}
		mj = mag(wid_j[XX], wid_j[YY], wid_j[ZZ]);
		mm = mag(wid_m[XX], wid_m[YY], wid_m[ZZ]);
		prod = dotp(wid_j[XX], wid_j[YY], wid_j[ZZ], wid_m[XX], wid_m[YY],
				wid_m[ZZ]) / (mm * mj);
		if (fabs(prod) < EPS) {
			*whperp = 1; /* width and height are perpend. to that on other fil*/
			return 1 == 1;
		} else if (fabs(fabs(prod) - 1) < EPS) {
			*whperp = 0;
			return 1 == 1;
		} else
			return 1 == 0;
	}
}

/* calculates direction of width if not specified */
void get_wid(FILAMENT *fil, double *wid) {

	double wx, wy, wz;
	double mag;

	if (fil->segm->widthdir != NULL) {
		wx = fil->segm->widthdir[XX];
		wy = fil->segm->widthdir[YY];
		wz = fil->segm->widthdir[ZZ];
	} else {
		/* default for width direction is in x-y plane perpendic to length*/
		/* so do cross product with unit z*/
		wx = -(fil->y[1] - fil->y[0]) * 1.0;
		wy = (fil->x[1] - fil->x[0]) * 1.0;
		wz = 0;
		if (fabs(wx / fil->segm->length) < EPS
				&& fabs(wy / fil->segm->length) < EPS) {
			/* if all of x-y is perpendic to length, then choose x direction */
			wx = 1.0;
			wy = 0;
		}
		mag = sqrt(wx * wx + wy * wy + wz * wz);
		wx = wx / mag;
		wy = wy / mag;
		wz = wz / mag;
	}
	wid[XX] = wx;
	wid[YY] = wy;
	wid[ZZ] = wz;
}

/* calculates direction of height */
void get_height(FILAMENT *fil, double *wid, double *height) {
	double wx = wid[XX];
	double wy = wid[YY];
	double wz = wid[ZZ];
	double hx, hy, hz, mag;

	hx = -wy * (fil->z[1] - fil->z[0]) + (fil->y[1] - fil->y[0]) * wz;
	hy = -wz * (fil->x[1] - fil->x[0]) + (fil->z[1] - fil->z[0]) * wx;
	hz = -wx * (fil->y[1] - fil->y[0]) + (fil->x[1] - fil->x[0]) * wy;
	mag = sqrt(hx * hx + hy * hy + hz * hz);
	hx = hx / mag;
	hy = hy / mag;
	hz = hz / mag;

	height[XX] = hx;
	height[YY] = hy;
	height[ZZ] = hz;
}

/* exact mutual inductance based on C. Hoer and C.Love, 
 Journal of the National Bureau of Standards-C,  Vol. 69C, p 127-137, 1965.
 */
double exact_mutual(FILAMENT *fil_j, FILAMENT *fil_m, int whperp, double *x_j,
		double *y_j, enum degen_type deg_j, enum degen_type deg_m)
/* x_j, y_j are unit vectors in the fil coord sys */
{
	double z_j[3]; /* unit vectors in the filament coord sys*/
	double origin[3];
	double dumb, ox, oy, oz, length;
	double a, b, c, d, l1, l2, l3, E, P, l3_1;
	double endx, endy, endz;
	int sign, sign2;
	double q[4], r[4], s[4];
	double totalM;
	int i, j, k;
	int a_deg, b_deg, c_deg, d_deg;
	extern int forced;

	/* the width direction of j will be the x direction in the filament
	 coordinate system */
	/*   get_wid(fil_j, x_j);  these are now passed as parms */

	/* the height direction will be the y direction */
	/*   get_height(fil_j, x_j, y_j);  these are now passed as parms */

	/*
	 z_j[XX] = fil_j->x[1] - fil_j->x[0];
	 z_j[YY] = fil_j->y[1] - fil_j->y[0];
	 z_j[ZZ] = fil_j->z[1] - fil_j->z[0];
	 dumb = mag(z_j[XX],z_j[YY],z_j[ZZ]);
	 */
	length = fil_j->length;
	z_j[XX] = fil_j->lenvect[XX] / length;
	z_j[YY] = fil_j->lenvect[YY] / length;
	z_j[ZZ] = fil_j->lenvect[ZZ] / length;

	a = fil_j->width;
	b = fil_j->height;
	if (whperp == 0) {
		c = fil_m->height;
		d = fil_m->width;
	} else {
		d = fil_m->height;
		c = fil_m->width;
	}

	ox = origin[XX] = fil_j->x[0] - x_j[XX] * a / 2 - y_j[XX] * b / 2;
	oy = origin[YY] = fil_j->y[0] - x_j[YY] * a / 2 - y_j[YY] * b / 2;
	oz = origin[ZZ] = fil_j->z[0] - x_j[ZZ] * a / 2 - y_j[ZZ] * b / 2;

	endx = fil_m->x[0] - ox;
	endy = fil_m->y[0] - oy;
	endz = fil_m->z[0] - oz;

	E = dotp(x_j[XX], x_j[YY], x_j[ZZ], endx, endy, endz) - d / 2;

	P = dotp(y_j[XX], y_j[YY], y_j[ZZ], endx, endy, endz) - c / 2;

	l3 = dotp(z_j[XX], z_j[YY], z_j[ZZ], endx, endy, endz);
	l3_1 = dotp(z_j[XX], z_j[YY], z_j[ZZ], fil_m->x[1] - ox, fil_m->y[1] - oy,
			fil_m->z[1] - oz);

	l1 = fil_j->segm->length;
	l2 = fil_m->segm->length;

	if (fabs(fabs(l3 - l3_1) - l2) / l2 > EPS) {
		fprintf(stderr, "Huh?  filament not expected length\n");
		exit(1);
	}

	if (l3 <= l3_1)
		sign = 1;
	else {
		sign = -1;
		l3 = l3_1;
	}

	a_deg = a / MAX(l1, b) < DEG_TOL;
	b_deg = b / MAX(l1, a) < DEG_TOL;
	c_deg = c / MAX(l2, d) < DEG_TOL;
	d_deg = d / MAX(l2, c) < DEG_TOL;

	if (forced) {
		if (deg_j == brick && deg_m == brick) {
			a_deg = d_deg = b_deg = c_deg = 0;
			fprintf(stderr, "fdbb");
		}
		if (deg_j == flat && deg_m == flat) {
			a_deg = d_deg = 0;
			b_deg = c_deg = 1;
			fprintf(stderr, "fdff ");
		}
		if (deg_j == flat && deg_m == skinny) {
			a_deg = c_deg = 0;
			b_deg = d_deg = 1;
			fprintf(stderr, "fdfs ");
		}
		if (deg_j == skinny && deg_m == flat) {
			a_deg = c_deg = 0;
			b_deg = d_deg = 1;
			fprintf(stderr, "fdsf ");
		}
		if (deg_j == skinny && deg_m == skinny) {
			a_deg = d_deg = 0;
			b_deg = c_deg = 1;
			fprintf(stderr, "fdss ");
		}
	}

	if (a_deg && b_deg && c_deg && d_deg) {
		/* two long filaments */
		totalM = fourfil(fil_j, fil_m) / (sign * MUOVER4PI);
	} else if (!a_deg && b_deg && c_deg && d_deg) {
		/* one flat and one long */
		totalM = tape_to_fil(E + d / 2, a, P - b / 2 + c / 2, l3, l1, l2) / a;
	} else if (!a_deg && b_deg && c_deg && !d_deg) {
		/* two flat */
		totalM = flat_to_flat_tape(E, a, d, P - b / 2 + c / 2, l3, l1, l2)
				/ (a * d);
	} else if (!a_deg && b_deg && !c_deg && d_deg) {
		/* one flat and one skinny */
		totalM = flat_to_skinny_tape(E + d / 2, a, P - b / 2, c, l3, l1, l2)
				/ (a * c);
	} else if (!a_deg && !b_deg && c_deg && d_deg) {
		totalM = brick_to_fil(E + d / 2, a, P + c / 2, b, l3, l1, l2) / (a * b);
	} else if (deg_j != brick || deg_m != brick) {
		fprintf(stderr, "Internal Error: Bad Degenerate filament a/l1<tol\n");
		fprintf(stderr, "Using fourfil() instead...\n");
		totalM = fourfil(fil_j, fil_m) / (sign * MUOVER4PI);
	} else
		/* all that are left are bricks */
		/* this is the full 3-D filament calculation, no degeneracies */
		totalM = brick_to_brick(E, a, d, P, b, c, l3, l1, l2) / (a * b * c * d);

	/*  moved to brick_to_brick()
	 fill_4(q, E,a,d);ads
	 fill_4(r, P,b,c);
	 fill_4(s, l3,l1,l2);

	 totalM = 0;

	 for(i = 0; i < 4; i++)
	 for(j = 0; j < 4; j++)
	 for(k = 0; k < 4; k++) {
	 sign2 = ( (i+j+k)%2 == 0 ? 1 : -1);
	 totalM += sign2*eval_eq(q[i],r[j],s[k], a);
	 }
	 */

#ifndef NO_ISNAN
#ifndef SGI
	if (isnan(totalM))
		fprintf(stderr, "Exact_mutual returned NaN! filaments %d, %d\n",
				fil_j->filnumber, fil_m->filnumber);
#endif
#endif

	return sign * MUOVER4PI * totalM;

	/*  return sign*MUOVER4PI*totalM/(a*b*c*d); */
}

void fill_4(double *vec, double E, double a, double d) {
	vec[0] = E - a;
	vec[1] = E + d - a;
	vec[2] = E + d;
	vec[3] = E;
}

double eval_eq(double x, double y, double z, double ref_len) {
	static double one_60 = 1.0 / 60.0;
	static double one_6 = 1.0 / 6.0;
	static double one_24 = 1.0 / 24.0;
	double retval;
	double len, xsq, ysq, zsq;
	int num_nearzero;
	int num_nearzero_sq;
	double one_over_ref_len;
	double one_over_ref_len_sq;

	one_over_ref_len = 1.0 / MAX(ref_len, (fabs(x) + fabs(y) + fabs(z)));
	one_over_ref_len_sq = SQUARE(one_over_ref_len);

	xsq = x * x;
	ysq = y * y;
	zsq = z * z;

	len = sqrt(xsq + ysq + zsq);

	retval = one_60 * len
			* (xsq * (xsq - 3 * ysq) + ysq * (ysq - 3 * zsq)
					+ zsq * (zsq - 3 * xsq));

	num_nearzero = nearzero(x * one_over_ref_len) + nearzero(y*one_over_ref_len)
	+ nearzero(z*one_over_ref_len);

	num_nearzero_sq = nearzero(
			xsq * one_over_ref_len_sq) + nearzero(ysq*one_over_ref_len_sq)
			+ nearzero(zsq*one_over_ref_len_sq);

	if (num_nearzero_sq < 2)
		retval += one_24
				* (log_term(x, xsq, ysq, zsq, len)
						+ log_term(y, ysq, xsq, zsq, len)
						+ log_term(z, zsq, xsq, ysq, len));

	if (num_nearzero < 1)
		retval -= one_6
				* (tan_term(x, y, z, zsq, len) + tan_term(x, z, y, ysq, len)
						+ tan_term(z, y, x, xsq, len));

	return retval;
}

double log_term(double x, double xsq, double ysq, double zsq, double len) {
	double retval;
	retval = ((6 * zsq - ysq) * ysq - zsq * zsq) * x
			* log((x + len) / sqrt(ysq + zsq));
	return retval;
}

double tan_term(double x, double y, double z, double zsq, double len) {
	double retval;
	retval = x * y * z * zsq * atan(x * y / (z * len));
	return retval;
}

/*  now a macro
 int nearzero(double x)
 {
 return (fabs(x) < EPS);
 }
 */

/* the following is code for the lookup table for common mutual terms */

/* The lookup table for gp's will be a global variable. Yuck! */

/* The lookup table head */
Table *table = NULL;
AllocInfo table_alloc;
AllocInfo double_alloc;

#define TABLE_SMALL OFF

/* lookup mutual term in table */

int lookup(FILAMENT *fil_j, FILAMENT *fil_m, int whperp, double *widj,
		double *heightj, double *retval, double *dims, int *dim_count,
		Table ***lastptr, int *p_num_dims)

/* widjm heightj arewidth and height vectors */
{
	int num_dims = NUM_DIMS;
	Table **s_table;

	if (whperp == 1) {
		*lastptr = NULL;
		return 0;
	}

#if TABLE_SMALL == ON

	/* this was an attempted shorter table but works worse than original! */

	/* for ground planes, look in global table */

	if (fil_j->segm->node[0]->gp != NULL && fil_m->segm->node[0] != NULL) {
		num_dims = *p_num_dims = 9;
		fill_dims(fil_j, fil_m, widj, heightj, dims,num_dims);
		s_table = &table; /* global table */
	}
	else if (fil_j->segm == fil_m->segm) {
		/* if fils on same seg, look up individual segment tables */
		num_dims = *p_num_dims = 6;
		fill_dims_seg(fil_j, fil_m, widj, heightj, dims, num_dims);
		s_table = &fil_j->segm->table;
	}
	else {
		/* otherwise don't look for it */
		*lastptr = NULL;
		return 0;
	}
#else
	/* this is the original lookup table which works better! */
	num_dims = *p_num_dims = 10;
	fill_dims(fil_j, fil_m, widj, heightj, dims, num_dims);
	s_table = &table; /* global table */

#endif

	if (find_dims(dims, num_dims, s_table, retval, dim_count, lastptr) == 0)
		return 0;
	else {
		if (dotprod(fil_j, fil_m) < 0)
			*retval *= -1;
		return 1;
	}

}

/* this fills the vector dims with the dimension information for fil_j
 and fil_m to later determine if the pair has been previous computed and
 is in the lookup table.  All dims should be positive!!
 */
void fill_dims(FILAMENT *fil_j, FILAMENT *fil_m, double *widthj,
		double *heightj, double *dims, int num_dims) {
	int j_first;
	int is_same;
	int i = 0;
	double x, y, z;
	double z_j[3], length;

#if TABLE_SMALL != ON
	if (fil_j->segm->node[0]->gp != NULL && fil_m->segm->node[0]->gp != NULL)
		dims[i++] = 1.0;
	else
		dims[i++] = 0.0;
#endif

	length = fil_j->length;
	/* vector along length */
	z_j[XX] = fil_j->lenvect[XX] / length;
	z_j[YY] = fil_j->lenvect[YY] / length;
	z_j[ZZ] = fil_j->lenvect[ZZ] / length;

	x = (fil_j->x[0] + fil_j->x[1]) - (fil_m->x[0] + fil_m->x[1]);
	y = (fil_j->y[0] + fil_j->y[1]) - (fil_m->y[0] + fil_m->y[1]);
	z = (fil_j->z[0] + fil_j->z[1]) - (fil_m->z[0] + fil_m->z[1]);

	/* compute center to center distances in fil_j coordinates */
	/* (putting height component first should sort the information by */
	/* plane (for multiple planes only)) */
	/*
	 dims[i++] = fabs(dotp(heightj[XX], heightj[YY], heightj[ZZ],x,y,z));
	 dims[i++] = fabs(dotp(widthj[XX], widthj[YY], widthj[ZZ],x,y,z));
	 dims[i++] = fabs(dotp(z_j[XX], z_j[YY], z_j[ZZ],x,y,z));
	 */

	if ((is_same = compare(fil_j->height, fil_m->height, EPS)) == -1)
		j_first = 1;
	else if (is_same == 1)
		j_first = 0;
	else if ((is_same = compare(fil_j->width, fil_m->width, EPS)) == -1)
		j_first = 1;
	else if (is_same == 1)
		j_first = 0;
	else if ((is_same = compare(fil_j->length, fil_m->length, EPS)) == -1)
		j_first = 1;
	else if (is_same == 1 || is_same == 0)
		j_first = 0;

	if (j_first == 1) {
		dims[i++] = fil_j->height;
		dims[i++] = fil_j->width;
		dims[i++] = fil_j->length;
		dims[i++] = fil_m->height;
		dims[i++] = fil_m->width;
		dims[i++] = fil_m->length;
	} else {
		dims[i++] = fil_m->height;
		dims[i++] = fil_m->width;
		dims[i++] = fil_m->length;
		dims[i++] = fil_j->height;
		dims[i++] = fil_j->width;
		dims[i++] = fil_j->length;
	}

	dims[i++] = fabs(dotp(heightj[XX], heightj[YY], heightj[ZZ], x, y, z));
	dims[i++] = fabs(dotp(widthj[XX], widthj[YY], widthj[ZZ], x, y, z));
	dims[i++] = fabs(dotp(z_j[XX], z_j[YY], z_j[ZZ], x, y, z));

	if (i != num_dims) {
		fprintf(stderr, "bad number for num_dims %d\n", num_dims);
		exit(1);
	}
}

/* this fills the vector dims with the dimension information for fil_j
 and fil_m to later determine if the pair has been previous computed and
 is in the lookup table.  It is for fils on the same seg only.
 */
void fill_dims_seg(FILAMENT *fil_j, FILAMENT *fil_m, double *widthj,
		double *heightj, double *dims, int num_dims) {
	int j_first;
	int is_same;
	int i = 0;
	double x, y, z;
	double z_j[3], length;

	x = (fil_j->x[0] + fil_j->x[1]) - (fil_m->x[0] + fil_m->x[1]);
	y = (fil_j->y[0] + fil_j->y[1]) - (fil_m->y[0] + fil_m->y[1]);
	z = (fil_j->z[0] + fil_j->z[1]) - (fil_m->z[0] + fil_m->z[1]);
	dims[i++] = fabs(dotp(widthj[XX], widthj[YY], widthj[ZZ], x, y, z));
	dims[i++] = fabs(dotp(heightj[XX], heightj[YY], heightj[ZZ], x, y, z));

	/* should always be zero for fils on same seg */
	/*  dims[i++] = fabs(dotp(z_j[XX], z_j[YY], z_j[ZZ],x,y,z)); */

	/* choose which fil to put first based on height, then width */
	if ((is_same = compare(fil_j->height, fil_m->height, EPS)) == -1)
		j_first = 1;
	else if (is_same == 1)
		j_first = 0;
	else if ((is_same = compare(fil_j->width, fil_m->width, EPS)) == -1)
		j_first = 1;
	else if (is_same == 1 || is_same == 0)
		j_first = 0;

	if (j_first == 1) {
		dims[i++] = fil_j->height;
		dims[i++] = fil_j->width;
		dims[i++] = fil_m->height;
		dims[i++] = fil_m->width;
	} else {
		dims[i++] = fil_m->height;
		dims[i++] = fil_m->width;
		dims[i++] = fil_j->height;
		dims[i++] = fil_j->width;
	}

	if (i != num_dims) {
		fprintf(stderr, "bad number for num_dims %d\n", num_dims);
		exit(1);
	}
}

/*
 int compare(double x, double y, double eps)
 {
 if ( (x==0 && y==0) || (fabs(x - y)/fabs(x+y) < eps))
 return 0;
 else if (x > y)
 return 1;
 else
 return -1;
 }
 */

int find_dims(double *dims, int num_dims, Table **a_table, double *retval,
		int *ret_dim_count, Table ***ret_lastptr) {
	Table *entry, **lastptr;
	int is_same, maybe_its_there;
	int dim_count, i;

	maybe_its_there = TRUE;
	dim_count = 0;
	/*  entry = table;  lastptr = &table;  the old code */
	entry = *a_table;
	lastptr = a_table;
	while (entry != NULL && maybe_its_there == TRUE) {
		if ((is_same = compare(dims[dim_count], entry->val, EPS)) == 1)
			/* not found */
			maybe_its_there = FALSE;
		else if (is_same == -1) {
			lastptr = &(entry->next_val);
			entry = entry->next_val;
		} else {
			if (dim_count < num_dims - 1) {
				dim_count++;
				lastptr = &(entry->u.next_dim);
				entry = entry->u.next_dim;
			} else {
				*retval = *(entry->u.mut_term);
				/*
				 printf("Found!    ");
				 for(i=0;i<num_dims;i++) printf("%6.3lg ",dims[i]);
				 printf("%13.6lg\n",*retval);
				 */
				return 1;
			}
		}
	}

	/* info for put_in_table to quickly insert in table */
	(*ret_lastptr) = lastptr;
	(*ret_dim_count) = dim_count;
	/*
	 printf("Not Found!");
	 for(i=0;i<num_dims;i++) printf("%6.3lg ",dims[i]);
	 printf("\n");
	 */
	return 0;
}

void put_in_table(FILAMENT *fil_j, FILAMENT *fil_m, int whperp, double mutterm,
		double *dims, int dim_count, Table **lastptr, int num_dims) {
	Table *entry;
	int i;

	if (lastptr == NULL)
		return;

	entry = (Table *) AllocAnEntry(&table_alloc);
	entry->next_val = (*lastptr);
	(*lastptr) = entry;
	entry->val = dims[dim_count++];
	lastptr = &(entry->u.next_dim);

	for (i = dim_count; i < num_dims; i++) {
		entry = (Table *) AllocAnEntry(&table_alloc);
		(*lastptr) = entry;
		entry->val = dims[i];
		entry->next_val = NULL;
		lastptr = &(entry->u.next_dim);
	}

	entry->u.mut_term = (double *) AllocAnEntry(&double_alloc);
	*(entry->u.mut_term) = fabs(mutterm);

	/*
	 printf("put in:   ");
	 for(i=0;i<num_dims;i++) printf("%6.3lg ",dims[i]);
	 printf("%13.6lg\n",fabs(mutterm));
	 */
}

/* no logic, just a nice number */
#define ALLOCBLOCK 256 

/* initialize info for allocating table so we can free it later */
void init_table(void) {
	table_alloc.size = sizeof(Table);
	table_alloc.blocksize = ALLOCBLOCK;
	table_alloc.elems_left = 0;
	table_alloc.head = NULL;
	double_alloc.size = sizeof(double);
	double_alloc.blocksize = ALLOCBLOCK;
	double_alloc.elems_left = 0;
	double_alloc.head = NULL;
}

int get_table_mem(void) {
	return MemoryForEntries(&table_alloc) + MemoryForEntries(&double_alloc);
}

void destroy_table(void) {
	DestroyEntries(&table_alloc);
	DestroyEntries(&double_alloc);
}

/* allocates table entries in blocks for more efficient memory */
char *AllocAnEntry(AllocInfo *allocptr) {
	int blocksize, size;
	AllocList *elem;

	if (allocptr->elems_left > 0) {
		allocptr->elems_left--;
		return (allocptr->next_elem += allocptr->size);
	} else {
		blocksize = allocptr->blocksize;
		elem = (AllocList *) malloc(sizeof(AllocList));
		elem->next = allocptr->head;
		allocptr->head = elem;
		elem->ptr = allocptr->next_elem = (char *) malloc(
				allocptr->size * blocksize);
		allocptr->elems_left = blocksize - 1;
		if (allocptr->next_elem == NULL) {
			fprintf(stderr,
					"Out of memory in AllocAnEntry. size = %d, block = %d\n",
					allocptr->size, blocksize);
			exit(1);
		}
		return allocptr->next_elem;
	}
}

void DestroyEntries(AllocInfo *allocinfo) {
	AllocList *lastelem, *head;

	head = allocinfo->head;

	while (head != NULL) {
		lastelem = head;
		free(head->ptr);
		head = head->next;
		free(lastelem);
	}

	allocinfo->head = NULL;
	allocinfo->elems_left = 0;
}

int MemoryForEntries(AllocInfo *allocptr) {
	AllocList *entry;
	int count = 0;

	for (entry = allocptr->head; entry != NULL; entry = entry->next)
		count++;

	return count * allocptr->blocksize * allocptr->size;
}

double brick_to_brick(double E, double a, double d, double P, double b,
		double c, double l3, double l1, double l2) {
	double q[4], r[4], s[4], totalM;
	int i, j, k, sign2;

	fill_4(q, E, a, d);
	fill_4(r, P, b, c);
	fill_4(s, l3, l1, l2);

	totalM = 0;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++) {
				sign2 = ((i + j + k) % 2 == 0 ? 1 : -1);
				totalM += sign2 * eval_eq(q[i], r[j], s[k], a);
			}

	return totalM;
}

double flat_to_flat_tape(double E, double a, double d, double P, double l3,
		double l1, double l2) {
	double q[4], s[4], totalM;
	int i, k, sign2;

	fill_4(q, E, a, d);
	fill_4(s, l3, l1, l2);

	totalM = 0;

	for (i = 0; i < 4; i++)
		for (k = 0; k < 4; k++) {
			sign2 = ((i + k) % 2 == 0 ? 1 : -1);
			totalM += sign2 * eval_eq_tape(q[i], P, s[k], a);
		}

	return totalM;
}

double eval_eq_tape(double x, double y, double z, double ref_len) {
	static double one_6 = 1.0 / 6.0;
	static double one_3 = 1.0 / 3.0;
	double retval;
	double len, xsq, ysq, zsq;
	double one_over_ref_len, one_over_ref_len_sq;
	int num_nearzero_sq;

	one_over_ref_len = 1.0 / MAX(ref_len, (fabs(x) + fabs(y) + fabs(z)));
	one_over_ref_len_sq = SQUARE(one_over_ref_len);

	xsq = x * x;
	ysq = y * y;
	zsq = z * z;

	len = sqrt(xsq + ysq + zsq);

	retval = -one_6 * len * (xsq - 2 * ysq + zsq);

	num_nearzero_sq = nearzero(
			xsq * one_over_ref_len_sq) + nearzero(ysq*one_over_ref_len_sq)
			+ nearzero(zsq*one_over_ref_len_sq);

	if (num_nearzero_sq < 2)
		retval += 0.5
				* ((xsq - ysq) * z * log(z + len)
						+ (zsq - ysq) * x * log(x + len));

	if (!nearzero(y * one_over_ref_len))
		retval -= x * y * z * atan(x * z / (y * len));

	return retval;
}

double flat_to_skinny_tape(double E, double a, double P, double c, double l3,
		double l1, double l2) {
	double q[2], r[2], s[4], totalM;
	int i, j, k, sign2;

	q[0] = E;
	q[1] = E - a;
	r[0] = P + c;
	r[1] = P;
	fill_4(s, l3, l1, l2);

	totalM = 0;

	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			for (k = 0; k < 4; k++) {
				sign2 = ((i + j + k) % 2 == 0 ? 1 : -1);
				totalM += sign2 * eval_eq_tape2(q[i], r[j], s[k], a);
			}

	return totalM;
}

double eval_eq_tape2(double x, double y, double z, double ref_len) {
	static double one_6 = 1.0 / 6.0;
	static double one_3 = 1.0 / 3.0;
	int num_nearzero;
	int num_nearzero_sq;
	double retval;
	double len, xsq, ysq, zsq;
	double one_over_ref_len, one_over_ref_len_sq;
	int nzxsq, nzysq, nzzsq;

	one_over_ref_len = 1.0 / MAX(ref_len, (fabs(x) + fabs(y) + fabs(z)));
	one_over_ref_len_sq = SQUARE(one_over_ref_len);

	xsq = x * x;
	ysq = y * y;
	zsq = z * z;

	len = sqrt(xsq + ysq + zsq);

	retval = -one_3 * len * x * y;

	nzxsq = nearzero(xsq * one_over_ref_len_sq);
	nzysq = nearzero(ysq * one_over_ref_len_sq);
	nzzsq = nearzero(zsq * one_over_ref_len_sq);

	if (!(nzzsq && nzysq))
		retval += (0.5 * zsq - one_6 * ysq) * y * log(x + len);
	if (!(nzzsq && nzxsq))
		retval += (0.5 * zsq - one_6 * xsq) * x * log(y + len);
	if (!(nzzsq || nzysq || nzxsq))
		retval += x * y * z * log(z + len);

	/*  retval += (0.5*zsq - one_6*ysq)*y*log(x + len)
	 + (0.5*zsq - one_6*xsq)*x*log(y + len)
	 + x*y*z*log(z + len);
	 */

	num_nearzero = nearzero(x * one_over_ref_len) + nearzero(y*one_over_ref_len)
	+ nearzero(z*one_over_ref_len);

	if (num_nearzero < 1)
		retval -= zsq * z * one_6 * tan_tape(x, y, z, len)
				+ 0.5 * z
						* (xsq * tan_tape(y, z, x, len)
								+ ysq * tan_tape(x, z, y, len));

	return retval;
}

double tan_tape(double x, double y, double z, double len) {
	return atan(x * y / (z * len));
}

double tape_to_fil(double E, double a, double P, double l3, double l1,
		double l2) {
	/* I have not implemented this degenerate case.  It should be done
	 by fourfil() and never get to this point */
	fprintf(stderr, "Hey, tape_to_fil should not have been called!\n");
	exit(1);
}

double brick_to_fil(double E, double a, double P, double b, double l3,
		double l1, double l2) {
	/* I have not implemented this degenerate case.  It should be done
	 by fourfil() and never get to this point */
	fprintf(stderr, "Hey, brick_to_fil should not have been called!\n");
	exit(1);
}

