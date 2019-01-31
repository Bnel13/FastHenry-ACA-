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


/* this file has some of the functions for mutual and self inductance
 calculation.  The rest are in joelself.c */

#include "induct.h"

/* these are missing in some math.h files */
#ifdef NO_ATANH
double atanh(double x) {return (0.5*log((1.0+x)/(1.0-x)));}
double asinh(double x) {return (log(x + sqrt(x*x+1)));}
#endif
#ifdef NO_ISNAN
int finite(double x) {return (1);}
#else
/* SRW -- finite() is deprecated */
#define finite(x) (!isnan(x) && !isinf(x))
/* #define finite isfinite */
#endif

/* SRW */
double mutual(FILAMENT*, FILAMENT*);
void print_infinity_warning(FILAMENT*, FILAMENT*);
void findfourfils(FILAMENT*, FILAMENT*);
double selfterm(FILAMENT*);
double mutualfil(FILAMENT*, FILAMENT*);
double magdiff2(FILAMENT*, int, FILAMENT*, int);
double mut_rect(double, double);
double dotprod(FILAMENT*, FILAMENT*);
double fourfil(FILAMENT*, FILAMENT*);
double parallel_fils(FILAMENT*, FILAMENT*, int, double*, double*, double);

int realcos_error = 0;

/* calculates mutual inductance of "filaments" with width and height */
/* as a combination of filament approximations                       */
/* some functions it uses are in dist_betw_fils.c */
double mutual(FILAMENT *fil_j, FILAMENT *fil_m) {
	double totalM;
	int i, j, ij;
	double aspect_j, aspect_m;
	static double cutoff = 3;
	int nfilsj, nfilsm;
	double dist, rj, rm, sum1, sum2;
	int parallel, whperp;
	int edge_par, num_dims;
	extern double **Gweight; /* gaussian quad weights. */
	/*  Filled in dist_betw_fils.c*/
	double widj[3], heightj[3];
	double dims[NUM_DIMS];
	Table **lastptr = NULL;//= NULL bap edit
	int dim_count;
	dist = dist_betw_fils(fil_j, fil_m, &parallel);
	rj = MAX(fil_j->width, fil_j->height) / 2.0;
	rm = MAX(fil_m->width, fil_m->height) / 2.0;

	if (MAX(rj,rm) * 100 < dist) {
		/* fils are far apart */
		/* if (dist != 0.0 && MAX(rj,rm)/dist < 0.2) printf("1"); */
		/* printf("\nmfil: %14.8le ", mutualfil(fil_j, fil_m)); */
		num_mutualfil++;

		totalM = mutualfil(fil_j, fil_m);

		if (!finite(totalM))
			print_infinity_warning(fil_j, fil_m);

		return totalM;
	} else {
		/*
		 aspect_j = aspectratio(fil_j);
		 aspect_m = aspectratio(fil_m);
		 */

		if (parallel == 1) {
			get_wid(fil_j, widj);
			get_height(fil_j, widj, heightj);
		} else if (fabs(vdotp(fil_j->lenvect, fil_m->lenvect))
				/ (fil_j->length * fil_m->length) < EPS) {
			num_perp++;
			/* fils are perpendicular */
			return 0.0;
		}

		edge_par = parallel == 1 && edges_parallel(fil_j, fil_m, widj, &whperp);
		if (edge_par)
			if (lookup(fil_j, fil_m, whperp, widj, heightj, &totalM, dims,
					&dim_count, &lastptr, &num_dims) == 1) {
				num_found++;
				return totalM;
			}
		if (edge_par && 2 * MAX(rj, rm) * 10 > dist) {
			/*(dist == 0.0 || MAX(rj,rm)/dist > 0.02*sqrt(MAX(aspect_j,aspect_m)))){*/

			/* fils are close enough to use exact integrals  (6/94) */
			totalM = parallel_fils(fil_j, fil_m, whperp, widj, heightj, dist);

			put_in_table(fil_j, fil_m, whperp, totalM, dims, dim_count, lastptr,
					num_dims);
			/* printf("ex: %14.8lg \n", totalM);*/

			if (!finite(totalM))
				print_infinity_warning(fil_j, fil_m);

			return totalM;
		}
		/* no longer: else if (aspect_j<cutoff&&aspect_m<cutoff || (rj+rm)<dist) */
		else {
			/* do 5 filament approximation to the filament */

			totalM = fourfil(fil_j, fil_m);

			if (edge_par)
				put_in_table(fil_j, fil_m, whperp, totalM, dims, dim_count,
						lastptr, num_dims);
			/* printf("4fil: %14.8lg\n ", totalM); */

			if (!finite(totalM))
				print_infinity_warning(fil_j, fil_m);

			return totalM;
		}
#if 1==0
		/* gaussian quadrature.  exact expression used instead */
		else {
			/* one of the fils has a high aspect ratio. Let's do 1-D Gaussian quad */
			nfilsj = MIN((int)aspect_j, MAXsubfils);
			nfilsm = MIN((int)aspect_m, MAXsubfils);
			findnfils(fil_j,subfilj, nfilsj);
			findnfils(fil_m,subfilm, nfilsm);

			sum1 = 0;
			for(i = 0; i < nfilsj; i++) {
				sum2 = 0;
				for(j = 0; j < nfilsm; j++)
				sum2 += Gweight[nfilsm][j]*mutualfil(&subfilj[i],&subfilm[j]);
				sum1 += Gweight[nfilsj][i]*sum2;
			}

			printf("gss: %14.8le ", sum1/4.0);
			return sum1/4.0; /* divide by 4.0 since Gquad is on [-1,1] */
		}
#endif
	}
}

void print_infinity_warning(FILAMENT *fil1, FILAMENT *fil2) {
	FILAMENT *fil;

	fprintf(stderr,
			"Severe warning: mutual inductance = infinity for two filaments:\n");

	fil = fil1;
	fprintf(stderr,
			"  fil 1: from %lg, %lg, %lg to %lg, %lg, %lg  width:%lg height: %lg\n",
			fil->x[0], fil->y[0], fil->z[0], fil->x[1], fil->y[1], fil->z[1],
			fil->width, fil->height);
	fil = fil2;
	fprintf(stderr,
			"  fil 2: from %lg, %lg, %lg to %lg, %lg, %lg  width:%lg height: %lg\n",
			fil->x[0], fil->y[0], fil->z[0], fil->x[1], fil->y[1], fil->z[1],
			fil->width, fil->height);
	fprintf(stderr,
			"Probably because there are overlapping but non-orthogonal segments in the input\n");
}

void findfourfils(FILAMENT *fil, FILAMENT *subfils) {
	double hx, hy, hz, mag, wx, wy, wz;
	int i;

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

	hx = -wy * (fil->z[1] - fil->z[0]) + (fil->y[1] - fil->y[0]) * wz;
	hy = -wz * (fil->x[1] - fil->x[0]) + (fil->z[1] - fil->z[0]) * wx;
	hz = -wx * (fil->y[1] - fil->y[0]) + (fil->x[1] - fil->x[0]) * wy;
	mag = sqrt(hx * hx + hy * hy + hz * hz);
	hx = hx / mag;
	hy = hy / mag;
	hz = hz / mag;

	/* all mutualfil needs are the filament coordinates and length */
	for (i = 0; i < 2; i++) {
		subfils[0].x[i] = fil->x[i] + fil->width * wx / 2;
		subfils[0].y[i] = fil->y[i] + fil->width * wy / 2;
		subfils[0].z[i] = fil->z[i] + fil->width * wz / 2;

		subfils[1].x[i] = fil->x[i] - fil->width * wx / 2;
		subfils[1].y[i] = fil->y[i] - fil->width * wy / 2;
		subfils[1].z[i] = fil->z[i] - fil->width * wz / 2;

		subfils[2].x[i] = fil->x[i] + fil->height * hx / 2;
		subfils[2].y[i] = fil->y[i] + fil->height * hy / 2;
		subfils[2].z[i] = fil->z[i] + fil->height * hz / 2;

		subfils[3].x[i] = fil->x[i] - fil->height * hx / 2;
		subfils[3].y[i] = fil->y[i] - fil->height * hy / 2;
		subfils[3].z[i] = fil->z[i] - fil->height * hz / 2;
	}

	for (i = 0; i < 4; i++)
		subfils[i].length = fil->length;
}

/* calculates selfinductance of rectangular filament */
/* it uses an exact expression for the 6-fold integral from Ruehli */
/* (the actual function comes is in file joelself.c */
double selfterm(FILAMENT *fil) {
	double approx, joelself;

	/*   approx = fil->length*MUOVER4PI*2
	 *(log(2*fil->length/(K*(fil->width + fil->height))) - 1); */
	joelself = MU0 * self(fil->width, fil->length, fil->height);
	/*   printf("Joel's function: %lg,  my approx: %lg\n",joelself, approx); */

	return joelself;
}

/* calculates the mutual inductance between two filaments */
/* from Grover, Chapter 7 */
double mutualfil(FILAMENT *fil1, FILAMENT *fil2) {
	double R, R1, R2, R3, R4, l, m;
	double cose, sine, u, v, d;
	double alpha, tmp1, tmp2, tmp3, tmp4, tmp5;
	double maxR, maxlength, sinsq, blah, minR;
	double scaleEPS, realcos;
	double vtemp;

	double omega, M;
	int signofM;

	double junk;

	/* for parallel filaments */
	double ux, uy, uz, Rx, Ry, Rz, x1_0, x1_1, x2_0, x2_1, magu;
	double dx, dy, dz, dotp, vx, vy, vz;
	double R1sq, R2sq, R3sq, R4sq, m2, l2, u2, v2, alpha2;

	R1sq = magdiff2(fil1, 1, fil2, 1);
	R2sq = magdiff2(fil1, 1, fil2, 0);
	R3sq = magdiff2(fil1, 0, fil2, 0);
	R4sq = magdiff2(fil1, 0, fil2, 1);
	R1 = sqrt(R1sq);
	R2 = sqrt(R2sq);
	R3 = sqrt(R3sq);
	R4 = sqrt(R4sq);

	maxR = minR = R1;
	if (R2 > maxR)
		maxR = R2;
	if (R3 > maxR)
		maxR = R3;
	if (R4 > maxR)
		maxR = R4;

	if (R2 < minR)
		minR = R2;
	if (R3 < minR)
		minR = R3;
	if (R4 < minR)
		minR = R4;

	l = fil1->length;
	m = fil2->length;
	maxlength = (l > m ? l : m);

	scaleEPS = minR / maxlength * 10;
	if (scaleEPS < 1)
		scaleEPS = 1;
	if (scaleEPS > 100)
		scaleEPS = 100;

	alpha = R4sq - R3sq + R2sq - R1sq;
	signofM = 1;

	/*
	 if (alpha < 0) {
	 signofM = -1;
	 tmp1 = R1;
	 R1 = R4;
	 R4 = tmp1;
	 tmp1 = R3;
	 R3 = R2;
	 R2 = tmp1;
	 alpha = -alpha;
	 }
	 */

	/* segments touching */
	if ((fabs(R1) < EPS) || (fabs(R2) < EPS) || (fabs(R3) < EPS)
			|| (fabs(R4) < EPS)) {
		if (fabs(R1) < EPS)
			R = R3;
		else if (fabs(R2) < EPS)
			R = R4;
		else if (fabs(R3) < EPS)
			R = R1;
		else
			R = R2;

		M = MUOVER4PI * 2 * (dotprod(fil1, fil2) / (l * m))
				* (l * atanh(m / (l + R)) + m * atanh(l / (m + R)));
		/* note: dotprod should take care of signofM */

		return M;
	}

	cose = alpha / (2 * l * m);
	if (fabs(cose) > 1)
		cose = (cose < 0 ? -1.0 : 1.0);
	blah = 1.0 - fabs(cose);

	/* let's use the real cosine */
	realcos = dotprod(fil1, fil2) / (l * m);
	/*realcos = dotprod(fil1, fil2)/(fil1->length*fil2->length);*/

	/* Segments are perpendicular! */
	if (fabs(realcos) < EPS)
		return 0.0;

	if (fabs((realcos - cose) / cose) > 0.1)
		if (realcos_error == 0) {
			fprintf(stderr, "Internal Warning: realcos = %lg,  cose = %lg\n",
					realcos, cose);
			fprintf(stderr,
					"  This may be due to two filaments that are separated \n\
by a distance 1e10 times their length\n");
			realcos_error = 1;
		}

	cose = realcos;

	/* filaments parallel */
	tmp1 = fabs(fabs(cose) - 1);
	/*  if ( fabs( fabs(cose) - 1) < scaleEPS*EPS*10.0) { */
	if (fabs(fabs(cose) - 1) < EPS) {
		/* determine a vector in the direction of d with length d */
		/* (d is the distance between the lines made by the filament */
		Rx = fil2->x[0] - fil1->x[0]; /* vector from fil1 to fil2 */
		Ry = fil2->y[0] - fil1->y[0];
		Rz = fil2->z[0] - fil1->z[0];
		ux = fil1->x[1] - fil1->x[0];
		uy = fil1->y[1] - fil1->y[0];
		uz = fil1->z[1] - fil1->z[0];
		magu = sqrt(ux * ux + uy * uy + uz * uz);
		ux = ux / magu; /* unit vector in direction of fil1 */
		uy = uy / magu;
		uz = uz / magu;

		dotp = ux * Rx + uy * Ry + uz * Rz; /* component of R in direction of fil1 */

		/* d vector is R vector without its component in the direction of fils */
		dx = Rx - dotp * ux;
		dy = Ry - dotp * uy;
		dz = Rz - dotp * uz;
		d = sqrt(dx * dx + dy * dy + dz * dz);

		/* let fil1 be the x axis, with node 0 being origin and u be */
		/* its positive direction */
		x1_0 = 0;
		x1_1 = l;

		/* x2_0 = dotprod( fil2.node0 - (fil1.node0 + d), u ) */
		/* (dotproduct just gives it correct sign) */
		vx = (fil2->x[0] - (fil1->x[0] + dx));
		vy = (fil2->y[0] - (fil1->y[0] + dy));
		vz = (fil2->z[0] - (fil1->z[0] + dz));
		x2_0 = vx * ux + vy * uy + vz * uz;
		vtemp = sqrt(vx * vx + vy * vy + vz * vz);

		/* same thing for x2_1 */
		vx = (fil2->x[1] - (fil1->x[0] + dx));
		vy = (fil2->y[1] - (fil1->y[0] + dy));
		vz = (fil2->z[1] - (fil1->z[0] + dz));
		x2_1 = vx * ux + vy * uy + vz * uz;

		if (fabs(
				(sqrt(vx * vx + vy * vy + vz * vz) - fabs(x2_1))
						/ (MAX(fabs(x2_0) + d, fabs(x2_1) + d))) > EPS) {
			printf("uh oh, segs don't seem parallel %lg\n",
					(sqrt(vx * vx + vy * vy * vz * vz) - fabs(x2_1)));
		}

		if (fabs(
				(vtemp - fabs(x2_0))
						/ (MAX(fabs(x2_0) + d, fabs(x2_1) + d))) > EPS) {
			printf("uh oh, segs don't seem parallel\n");
		}

		/*
		 if (x2_1 < x2_0) {
		 tmp1 = x2_0;
		 x2_0 = x2_1;
		 x2_1 = tmp1;
		 if (signofM != -1)
		 printf("uh oh, inconsistent direction fil1x = %lf %lf fil2x=%lf %lf\n"
		 ,x1_0, x1_1, x2_0, x2_1);
		 }
		 */

		if (fabs(d) < EPS) { /* collinear! */
			/* SRW -- For whatever reason, the gcc in Red Hat Linux 6.0 can't
			 * evaluate the expression below.  It returns a NaN, which completely
			 * screws up the entire run.  This is a compiler bug.  Breaking it
			 * up a bit seems to have fixed the problem.
			 */
			/*
			 M = MUOVER4PI*(fabs(x2_1 - x1_0)*log(fabs(x2_1 - x1_0))
			 - fabs(x2_1 - x1_1)*log(fabs(x2_1 - x1_1))
			 - fabs(x2_0 - x1_0)*log(fabs(x2_0 - x1_0))
			 + fabs(x2_0 - x1_1)*log(fabs(x2_0 - x1_1)) );
			 return M;
			 */
			M = fabs(x2_1 - x1_0) * log(fabs(x2_1 - x1_0));
			M -= fabs(x2_1 - x1_1) * log(fabs(x2_1 - x1_1));
			M -= fabs(x2_0 - x1_0) * log(fabs(x2_0 - x1_0));
			M += fabs(x2_0 - x1_1) * log(fabs(x2_0 - x1_1));
			return (MUOVER4PI * M);

		} /* end collinear */

		M = MUOVER4PI
				* (mut_rect(x2_1 - x1_1, d) - mut_rect(x2_1 - x1_0, d)
						- mut_rect(x2_0 - x1_1, d) + mut_rect(x2_0 - x1_0, d));

		return M;
	} /* end if parallel filaments */

	/* the rest if for arbitrary filaments */

	l2 = l * l;
	m2 = m * m;
	alpha2 = alpha * alpha;

	u = l * (2 * m2 * (R2sq - R3sq - l2) + alpha * (R4sq - R3sq - m2))
			/ (4 * l2 * m2 - alpha2);
	v = m * (2 * l2 * (R4sq - R3sq - m2) + alpha * (R2sq - R3sq - l2))
			/ (4 * l2 * m2 - alpha2);

	u2 = u * u;
	v2 = v * v;

	d = (R3sq - u2 - v2 + 2 * u * v * cose);
	if (fabs(
			d / (R3sq + u2 + v2 + 1)
					* (maxlength * maxlength / (maxR * maxR))) < EPS)
		d = 0.0;
	d = sqrt(d);

	sinsq = 1.0 - cose * cose;
	if (fabs(sinsq) < EPS)
		sinsq = 0.0;
	sine = sqrt(sinsq);
	tmp1 = d * d * cose;
	tmp2 = d * sine;
	tmp3 = sine * sine;

	if (fabs(d) < EPS)
		omega = 0.0; /* d is zero, so it doesn't matter */
	else
		omega = atan2((tmp1 + (u + l) * (v + m) * tmp3), (tmp2 * R1))
				- atan2((tmp1 + (u + l) * v * tmp3), (tmp2 * R2))
				+ atan2((tmp1 + u * v * tmp3), (tmp2 * R3))
				- atan2((tmp1 + u * (v + m) * tmp3), (tmp2 * R4));

	tmp4 = ((u + l) * atanh(m / (R1 + R2)) + (v + m) * atanh(l / (R1 + R4))
			- u * atanh(m / (R3 + R4)) - v * atanh(l / (R2 + R3)));

	if (fabs(sine) < 1e-150)
		tmp5 = 0.0;
	else
		tmp5 = omega * d / sine;

	M = MUOVER4PI * cose * (2 * tmp4 - tmp5);

	return M;
}

/*
 double magdiff(FULAMENT *fil1, int node1, FULAMENT *fil2, int node2)
 {
 return sqrt( SQUARE(fil1->x[node1] - fil2->x[node2])
 +SQUARE(fil1->y[node1] - fil2->y[node2])
 +SQUARE(fil1->z[node1] - fil2->z[node2])
 );
 }
 */

double magdiff2(FILAMENT *fil1, int node1, FILAMENT *fil2, int node2) {
	return ( SQUARE(fil1->x[node1] - fil2->x[node2])
			+ SQUARE(fil1->y[node1] - fil2->y[node2])
			+ SQUARE(fil1->z[node1] - fil2->z[node2]));
}

/* this gives the mutual inductance of two filaments who represent */
/* opposite sides of a rectangle */

double mut_rect(double len, double d) {
	double temp, temp1;

	temp = sqrt(len * len + d * d);
	temp1 = len * asinh(len / d);
	return temp - temp1;
}

/* returns the dotproduct of the vector from node0 to node1 of fil1 */
/* with that of fil2 */
double dotprod(FILAMENT *fil1, FILAMENT *fil2) {
	return ((fil1->x[1] - fil1->x[0]) * (fil2->x[1] - fil2->x[0])
			+ (fil1->y[1] - fil1->y[0]) * (fil2->y[1] - fil2->y[0])
			+ (fil1->z[1] - fil1->z[0]) * (fil2->z[1] - fil2->z[0]));
}

double fourfil(FILAMENT *fil_j, FILAMENT *fil_m) {
	FILAMENT subfilj[MAXsubfils], subfilm[MAXsubfils];
	double totalM;
	int i;

	/* approximate 'filament' with width and length as a combination */
	/* of four filaments on the midpoints of the edges               */
	/* Known as the Rayleigh Quadrature formula. Grover p.11         */
	findfourfils(fil_j, subfilj);
	findfourfils(fil_m, subfilm);
	totalM = 0.0;
	for (i = 0; i < 4; i++)
		totalM += mutualfil(fil_j, &subfilm[i]);
	for (i = 0; i < 4; i++)
		totalM += mutualfil(fil_m, &subfilj[i]);
	totalM += -2.0 * mutualfil(fil_j, fil_m);

	totalM = totalM / 6.0;

	/* printf("4: %14.8le ",totalM); */
	num_fourfil++;

	return totalM;
}

double parallel_fils(FILAMENT *fil_j, FILAMENT *fil_m, int whperp, double *x_j,
		double *y_j, double dist)
/* double *x_j, *y_j;  unit vectors in the fil coord sys */
{
	enum degen_type deg_j, deg_m;

	/* find degenerate dimensions */
	deg_j = find_deg_dims(fil_j);
	deg_m = find_deg_dims(fil_m);

	if (deg_j == brick && deg_m == brick) {
		/* no degenerate dimensions, both are bricks */
		num_exact_mutual++;

		return exact_mutual(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m);

		/* fprintf(stderr,"Nondegenerate: fil %d to fil %d: %13.6lg\n",
		 fil_j->filnumber,
		 fil_m->filnumber,
		 exact_mutual(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m)); */
	} else {

		return compute_for_degenerate(fil_j, fil_m, whperp, x_j, y_j, deg_j,
				deg_m, dist);

		/* fprintf(stderr,"  degenerate: fil %d to fil %d: %13.6lg\n",
		 fil_j->filnumber,
		 fil_m->filnumber,
		 compute_for_degenerate(fil_j, fil_m, whperp, x_j, y_j,
		 deg_j, deg_m, dist));
		 */

	}

}
