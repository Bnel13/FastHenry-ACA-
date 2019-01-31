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


/* The following are support functions for holes in a ground plane. */
/*  For user defined hole functions, alter hole_user1() - hole_user10() */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "induct.h"

/* SRW */
int is_next_word(char*, char*);
int is_hole(NODES*);
HoleList *make_holelist(HoleList*, char*, double, double, double, double, int*);
int skipspace(char*);
int eos(char);
void hole_error(char*, char*, HoleList*);
int is_one_of(char, char*);
void delete_node(NODES*);
void make_holes(HoleList*, GROUNDPLANE*);
void hole_point(HoleList*, GROUNDPLANE*, double, double, double, double);
void hole_rect(HoleList*, GROUNDPLANE*, double, double, double, double);
void hole_circle(HoleList*, GROUNDPLANE*, double, double, double, double);
void hole_user1(HoleList*, GROUNDPLANE*, double, double, double, double);
void hole_user2(HoleList*, GROUNDPLANE*, double, double, double, double);
void hole_user3(HoleList*, GROUNDPLANE*, double, double, double, double);
void hole_user4(HoleList*, GROUNDPLANE*, double, double, double, double);
void hole_user5(HoleList*, GROUNDPLANE*, double, double, double, double);
void hole_user6(HoleList*, GROUNDPLANE*, double, double, double, double);
void hole_user7(HoleList*, GROUNDPLANE*, double, double, double, double);

/* returns TRUE if the next word in line matches str. */
/* (strips leading spaces from line, but not str) */
int is_next_word(char *str, char *line) {

	/* skip leading spaces */
	while (isspace(*line) && *line != '\0')
		line++;

	if (*line == '\0') /* anything left? */
		return FALSE;
	else {
		while (*str != '\0' && *str == *line) {
			str++;
			line++;
		}
		if (*str == '\0')
			return TRUE;
		else
			return FALSE;
	}
}

/* has the node been removed because it is part of a ground plane hole? */
int is_hole(NODES *node) {
	return (node->type == GPHOLE);
}

#define MAXNAME 20

/* this turns the info in line into a HoleList element */
HoleList *make_holelist(HoleList *hole_head, char *line, double units,
		double relx, double rely, double relz, int *skip) {
	char name[MAXNAME];
	int skip1, skip2, skip3, skip4;
	HoleList *holep;
	char *linep;
	double *vals;

	holep = (HoleList *) Gmalloc(sizeof(HoleList));

	name[MAXNAME - 1] = '\0';
	sscanf(line, "%19s%n", name, &skip1);
	line += skip1;

	if (strcmp("hole", name) != 0) {
		fprintf(stderr, "Internal Error:  %s is not the word 'hole'\n", name);
		exit(1);
	}

	name[MAXNAME - 1] = '\0';
	sscanf(line, "%19s%n", name, &skip2);
	line += skip2;

	if (name[MAXNAME - 1] != '\0')
		printf("Warning: hole function '%s' truncated to 19 chars\n", name);

	holep->func = (char *) Gmalloc((strlen(name) + 1) * sizeof(char));
	strcpy(holep->func, name);
	holep->units = units;
	holep->relx = relx;
	holep->rely = rely;
	holep->relz = relz;

	skip3 = 0;
	while (isspace(*line) && *line != '\0') {
		line++;
		skip3++;
	}

	if (*line != '(')
		hole_error("Values for hole must start with '('", line, holep);

	/* let's count how many values we have first */
	linep = line;
	holep->numvals = 0;

	while (*linep != ')' && !eos(*linep)) {
		linep++;
		while (isspace(*linep) || is_one_of(*linep, "1234567890.e+-"))
			linep++;
		holep->numvals++;
	}
	if (*linep != ')')
		hole_error("Hole values did not end with ')'", line, holep);

	holep->vals = (double *) Gmalloc(holep->numvals * sizeof(double));

	skip3 += linep - line + 1;

	/* now let's read them in */
	linep = line + 1; /* skip the ( */
	vals = holep->vals;

	while (*linep != ')') {
		if (sscanf(linep, "%lf%n", vals, &skip4) != 1)
			hole_error("Couldn't read value starting at:", linep, holep);
		else {
			linep += skip4;
			vals++;
			linep += skipspace(linep);
			if (*linep == ',')
				linep++;
		}
	}

	*skip = skip1 + skip2 + skip3;

	holep->next = hole_head;
	return holep;

}

int skipspace(char *line) {
	int skip = 0;

	while (isspace(*line)) {
		skip++;
		line++;
	}

	return skip;
}

/* End of string test */
int eos(char chr) {
	return (chr == '\0');
}

void hole_error(char *errstr, char *line, HoleList *holep) {
	fprintf(stderr, "While reading in a %s hole:\n", holep->func);
	fprintf(stderr, "%s\n%s", errstr, line);
	fprintf(stderr, "The hole format is 'hole <type> (val1,val2,...)'\n");
	fprintf(stderr,
			"Where type is a string, and valn is a floating point value\n");
	exit(1);
}

int is_one_of(char letter, char *one_of) {
	for (; *one_of != '\0'; one_of++)
		if (letter == *one_of)
			return (1 == 1);

	return 1 == 0;
}

void delete_node(NODES *node) {
	node->type = GPHOLE;
}

/* This calls the functions which make the holes */
void make_holes(HoleList *holep, GROUNDPLANE *gp) {
	double relx = holep->relx;
	double rely = holep->rely;
	double relz = holep->relz;
	double units = holep->units;

	if (strncmp("rect", holep->func, 4) == 0)
		hole_rect(holep, gp, relx, rely, relz, units);
	else if (strncmp("circle", holep->func, 6) == 0)
		hole_circle(holep, gp, relx, rely, relz, units);
	else if (strncmp("point", holep->func, 5) == 0)
		hole_point(holep, gp, relx, rely, relz, units);
	else if (strncmp("user1", holep->func, 5) == 0)
		hole_user1(holep, gp, relx, rely, relz, units);
	else if (strncmp("user2", holep->func, 5) == 0)
		hole_user2(holep, gp, relx, rely, relz, units);
	else if (strncmp("user3", holep->func, 5) == 0)
		hole_user3(holep, gp, relx, rely, relz, units);
	else if (strncmp("user4", holep->func, 5) == 0)
		hole_user4(holep, gp, relx, rely, relz, units);
	else if (strncmp("user5", holep->func, 5) == 0)
		hole_user5(holep, gp, relx, rely, relz, units);
	else if (strncmp("user6", holep->func, 5) == 0)
		hole_user6(holep, gp, relx, rely, relz, units);
	else if (strncmp("user7", holep->func, 5) == 0)
		hole_user7(holep, gp, relx, rely, relz, units);
	else
		hole_error("Unknown type of hole", "", holep);
}

/* This makes a hole by removing one node nearest to the (x,y,z) given
 by holes->vals;
 */
void hole_point(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
	double *vals = holep->vals;
	int i1, j1;

	if (holep->numvals != 3)
		hole_error("Exactly 3 values required for a point hole.", "", holep);

	delete_node(
			find_nearest_gpnode(vals[0] * units + relx, vals[1] * units + rely,
					vals[2] * units + relz, gp, &i1, &j1));
}

/* The following function creates a hole which is in the shape of 
 a rectangle whose edges are parallel to the ground plane edges.
 The 6 values in holep->vals are two (x,y,z) pairs which represent
 opposite corners of the rectangle.
 */
void hole_rect(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
	NODES *node1, *node2;
	int i1, j1, i2, j2;
	double *vals = holep->vals;
	int i_beg, i_end, j_beg, j_end, i, j;

	if (holep->numvals != 6)
		hole_error("Exactly 6 values required for a square hole.", "", holep);

	node1 = find_nearest_gpnode(vals[0] * units + relx, vals[1] * units + rely,
			vals[2] * units + relz, gp, &i1, &j1);

	node2 = find_nearest_gpnode(vals[3] * units + relx, vals[4] * units + rely,
			vals[5] * units + relz, gp, &i2, &j2);

	if (i1 <= i2) {
		i_beg = i1;
		i_end = i2;
	} else {
		i_beg = i2;
		i_end = i1;
	}

	if (j1 <= j2) {
		j_beg = j1;
		j_end = j2;
	} else {
		j_beg = j2;
		j_end = j1;
	}

	for (i = i_beg; i <= i_end; i++)
		for (j = j_beg; j <= j_end; j++)
			delete_node(gp->pnodes[i][j]);

}

/* The following function makes a hole in the shape of a circle.
 It takes 4 values in holep->vals.  The first three are (x,y,z) of
 the center and the last is the radius, R.
 */
void hole_circle(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
	NODES *node1, *node2, *nodec;
	int i1, j1, i2, j2;
	double center[3];
	int i_beg, i_end, j_beg, j_end, i, j, ic, jc;
	double R;
	double ref1, ref2;
	double edge2;
	int k_up, k_down, p_left, p_right;

	if (holep->numvals != 4)
		hole_error("Exactly 4 values required for a circular hole.",
				"(x,y,z,radius)", holep);

	center[0] = holep->vals[0] * units + relx;
	center[1] = holep->vals[1] * units + rely;
	center[2] = holep->vals[2] * units + relz;

	R = holep->vals[3] * units;

	/* find node nearest to center as a place to start (a 'reference')*/
	nodec = find_nearest_gpnode(center[0], center[1], center[2], gp, &ic, &jc);

	/* The plane has it's own coordinate system with unit vectors (ux1,uy1,uz1)
	 and (ux2, uy2, uz2) which are directed along the edges of the plane.
	 Let's do everything in terms of them */

	/*find component in the u1 and u2 direction of vector from center to nodec */
	ref1 = dotp(nodec->x - center[0], nodec->y - center[1],
			nodec->z - center[2], gp->ux1, gp->uy1, gp->uz1);
	ref2 = dotp(nodec->x - center[0], nodec->y - center[1],
			nodec->z - center[2], gp->ux2, gp->uy2, gp->uz2);

	if (ref1 * ref1 + ref2 * ref2 > R) {
		fprintf(stderr,
				"Warning: Discretization too coarse for circle of radius %lf.\n No circlular hole made.\n",
				holep->vals[3]);
		return;
	}

	/* We will remove nodes inside the circle in strips running in u2 direction*/
	k_up = (R - ref1) / gp->d1;
	k_down = (-R - ref1) / gp->d1;

	if (k_up < 0 || k_down > 0) {
		printf("Internal Error: k_up or k_down of the wrong sign\n");
		exit(1);
	}

	/* loop through all the u2 directed strips */
	for (i = k_down; i <= k_up; i++) {
		/* find the u2 component of a point on the strip at u1 = ref1 + i*d1 */
		edge2 = sqrt(R * R - SQUARE(ref1 + i * gp->d1));
		/* count number of cells to the right and left */
		p_right = (edge2 - ref2) / gp->d2;
		p_left = (-edge2 - ref2) / gp->d2;

		if (-edge2 - ref2 < 0 && edge2 - ref2 > 0)
			for (j = p_left; j <= p_right; j++)
				delete_node(gp->pnodes[i + ic][j + jc]);
	}
}

void hole_user1(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
}

void hole_user2(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
}

void hole_user3(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
}

void hole_user4(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
}

void hole_user5(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
}

void hole_user6(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
}

void hole_user7(HoleList *holep, GROUNDPLANE *gp, double relx, double rely,
		double relz, double units) {
}

