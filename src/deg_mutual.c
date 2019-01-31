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



/* This file contains functions for computing the mutual inductance 
 between fils which have one or more degenerate dimensions. For example,
 it's width = 10^4 * height */

#include "induct.h"

#define LEN 4
#define WID 2
#define HEIGHT 1

/* SRW */
enum degen_type find_deg_dims(FILAMENT*);
double compute_for_degenerate(FILAMENT*, FILAMENT*, int, double*, double*,
		enum degen_type, enum degen_type, double);
void setup_tape_to_tape(FILAMENT*, FILAMENT*, int, double*, double*,
		enum degen_type, enum degen_type, FILAMENT*, FILAMENT*, double**,
		double**);
double do_tape_to_brick(FILAMENT*, FILAMENT*, int, double*, double*,
		enum degen_type, enum degen_type);

enum degen_type find_deg_dims(FILAMENT *fil) {
	double max;

	max = MAX(fil->length, fil->width);
	max = MAX(max, fil->height);

	return (fil->length / max < DEG_TOL) * LEN
			+ (fil->width / max < DEG_TOL) * WID
			+ (fil->height / max < DEG_TOL) * HEIGHT;
}

double compute_for_degenerate(FILAMENT *fil_j, FILAMENT *fil_m, int whperp,
		double *x_j, double *y_j, enum degen_type deg_j, enum degen_type deg_m,
		double dist)
/* double *x_j, *y_j;  unit vectors in the fil coord sys */
{

	FILAMENT nfil_j, nfil_m; /* new temp fils */
	double *nx_j, *ny_j;

	if (deg_j == brick && deg_m == brick) {
		/* neither is degenerate, this shouldn't happen */
		fprintf(stderr, "Hey, compute_degenerate was called, impossible!\n");
		exit(1);
	}

	if ((deg_j == flat || deg_j == skinny)
			&& (deg_m == flat || deg_m == skinny)) {
		setup_tape_to_tape(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m,
				&nfil_j, &nfil_m, &nx_j, &ny_j);
		return exact_mutual(&nfil_j, &nfil_m, whperp, nx_j, ny_j, deg_j, deg_m);
	} else if ((deg_m == brick && (deg_j == flat || deg_j == skinny))
			|| (deg_j == brick && (deg_m == flat || deg_m == skinny)))
		return do_tape_to_brick(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m);
	else if (deg_j == too_long && deg_m == too_long)
		return fourfil(fil_j, fil_m);
	else if (deg_j == too_long || deg_j == too_long)
		return fourfil(fil_j, fil_m);
	else
		return fourfil(fil_j, fil_m);

}

void setup_tape_to_tape(FILAMENT *fil_j, FILAMENT *fil_m, int whperp,
		double *x_j, double *y_j, enum degen_type deg_j, enum degen_type deg_m,
		FILAMENT *nfil_j, FILAMENT *nfil_m, double **nx_j, double **ny_j) {

	if (deg_j == flat) {
		*nfil_j = *fil_j;
		*nfil_m = *fil_m;
		*nx_j = x_j;
		*ny_j = y_j;
	} else if (deg_j == skinny) {
		/* turn skinny into flat orientation */
		*nfil_j = *fil_j;
		*nfil_m = *fil_m;
		/* swap coord sys */
		*ny_j = x_j;
		*nx_j = y_j;
		/* swap height and width */
		nfil_j->width = fil_j->height;
		nfil_j->height = fil_j->width;
		nfil_m->width = fil_m->height;
		nfil_m->height = fil_m->width;
	}
}

double do_tape_to_brick(FILAMENT *fil_j, FILAMENT *fil_m, int whperp,
		double *x_j, double *y_j, enum degen_type deg_j, enum degen_type deg_m) {

	FILAMENT nfil_j, nfil_m;
	double *nx_j = NULL, *ny_j = NULL, *dR;
	double wid_brick[3], hei_brick[3], orig_x[2], orig_y[2], orig_z[2];
	double x_flat[3], y_flat[3];
	double small_dim, sum;
	int i, j, gpoints;
	extern double **Gweight, **Gpoint; /* gaussian quad weights. */
	enum degen_type ndeg_j, ndeg_m;

	/*
	 if ( deg_m == brick && (deg_j == flat || deg_j == skinny)
	 || deg_j == brick && (deg_m == flat || deg_m == skinny))
	 return do_tape_to_brick(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m);
	 */

	if (deg_j == flat) {
		nfil_j = *fil_j;
		nfil_m = *fil_m;
		nx_j = x_j;
		ny_j = y_j;
		get_wid(fil_m, wid_brick);
		get_height(fil_m, wid_brick, hei_brick);
	} else if (deg_j == skinny) {
		/* turn skinny into flat orientation */
		nfil_j = *fil_j;
		nfil_m = *fil_m;
		/* swap coord sys */
		ny_j = x_j;
		nx_j = y_j;
		/* swap height and width */
		nfil_j.width = fil_j->height;
		nfil_j.height = fil_j->width;
		nfil_m.width = fil_m->height;
		nfil_m.height = fil_m->width;
		/* get them swapped */
		get_wid(fil_m, hei_brick);
		get_height(fil_m, hei_brick, wid_brick);
	} else if (deg_j == brick) {
		/* swap j and m */
		nfil_j = *fil_m;
		nfil_m = *fil_j;
		get_wid(fil_m, x_flat);
		get_height(fil_m, x_flat, y_flat);

		if (deg_m == flat) {
			nx_j = x_flat;
			ny_j = y_flat;
			for (i = 0; i < 3; i++) {
				wid_brick[i] = x_j[i];
				hei_brick[i] = y_j[i];
			}

		} else {
			nx_j = y_flat;
			ny_j = x_flat;
			nfil_j.width = fil_m->height;
			nfil_j.height = fil_m->width;
			nfil_m.width = fil_j->height;
			nfil_m.height = fil_j->width;
			for (i = 0; i < 3; i++) {
				wid_brick[i] = y_j[i];
				hei_brick[i] = x_j[i];
			}

		}
	}

	/* store original brick position */
	for (i = 0; i < 2; i++) {
		orig_x[i] = nfil_m.x[i];
		orig_y[i] = nfil_m.y[i];
		orig_z[i] = nfil_m.z[i];
	}

	if (nfil_m.width > nfil_m.height) {
		/* the height direction will be done discretely */
		small_dim = nfil_m.height / 2;
		nfil_m.height = 0;
		dR = hei_brick;
		if (whperp == 0) /* useful for testing only. if forced == 1 */
			ndeg_m = flat;
		else
			ndeg_m = skinny;
	} else {
		/* the width direction will be done discretely */
		small_dim = nfil_m.width / 2;
		nfil_m.width = 0;
		dR = wid_brick;
		if (whperp == 0)
			ndeg_m = skinny;
		else
			ndeg_m = flat;
	}

	/* if forced == 1, then setting ndeg_j matters */
	ndeg_j = flat;
	nfil_j.height = 0.0; /* insure we use the middle of filament x-section*/

	gpoints = 3;
	/* now do gaussian quadrature of tape_to_tape to approximate */
	sum = 0;
	for (i = 0; i < gpoints; i++) {
		for (j = 0; j < 2; j++) {
			nfil_m.x[j] = orig_x[j] + dR[XX] * small_dim * Gpoint[gpoints][i];
			nfil_m.y[j] = orig_y[j] + dR[YY] * small_dim * Gpoint[gpoints][i];
			nfil_m.z[j] = orig_z[j] + dR[ZZ] * small_dim * Gpoint[gpoints][i];
		}
		sum += Gweight[gpoints][i]
				* exact_mutual(&nfil_j, &nfil_m, whperp, nx_j, ny_j, ndeg_j,
						ndeg_m);
	}

	return sum / 2.0;
}

