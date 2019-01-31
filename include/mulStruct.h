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



/* # ***** sort to /src/header
 # ***** */
#include "patran.h"		/* for neutral file interface */

struct surface { /* a surface file and its permittivities */
	int type; /* CONDTR, DIELEC or BOTH */
	double trans[3]; /* translation vector applied to surface */
	double ref[3]; /* reference point for normal alignment */
	int ref_inside; /* TRUE=>ref point inside (on - side of) surf*/
	int end_of_chain; /* TRUE=>end o/surf chain all w/same cond#'s */
	char *title; /* title inside surface file */
	char *name; /* surface file name (neutral or quick fmat) */
	char *group_name; /* name of cond group (if any) surf part of */
	double outer_perm; /* relative permitivity on side n pnts into
	 (n is taken as the normal to 1st panel) */
	double inner_perm; /* relative permitivity on other side
	 (always zero for type = CONDTR) */
	struct charge *panels; /* linked list of panels on this surface */
	int num_panels; /* number of panels (incl dummies) this surf */
	int num_dummies; /* number of dummy panels on this surf */
	struct surface *next; /* linked list pointers */
	struct surface *prev;
};
typedef struct surface surface;

struct charge { /* point charge */
	struct charge *next; /* Next charge in linked list. */
	double corner[4][3]; /* Corner point coordinates. */
	int shape; /* 4=quad panel, 3 = triangular panel. */
	int index; /* charge value index in q array */
	double X[3], Y[3], Z[3]; /* Def coord system, Z is normal direction. */
	double max_diag; /* Longest diagonal of panel. */
	double min_diag; /* Shortest diagonal. */
	double length[4]; /* Edge lengths. */
	double area; /* Area of two triangluar regions. */
	double x, y, z; /* Centroid of the quadlrilateral.  */
	double moments[16]; /* Moments of the panel. */
	double *multipole; /* Temporarily Stores Q2M. */
	int cond; /* Conductor number. */
	int order; /* Multipole order. */
	/* things needed to do electric field evaluation by divided difference */
	int dummy; /* TRUE => dummy panel for elec field eval */
	struct surface *surf; /* surface that contains this panel */
	struct charge *pos_dummy; /* eval pnt w/pos displacement from x,y,z */
	struct charge *neg_dummy; /* eval pnt w/neg displacement from x,y,z */

	struct Filament *fil; /* real filament that this corresponds to */
};
typedef struct charge charge;

struct cube {
	/* Definition variables. */
	int index; /* unique index */
	int level; /* 0 => root */
	double x, y, z; /* Position of cube center. */
	int j, k, l; /* cube is cubes[level][j][k][l] */
	int flag; /* used for marking for tree walks */

	/* Upward Pass variables. */
	int mul_exact; /* TRUE => do not build a multipole expansn */
	struct cube *mnext; /* Ptr to next cube on which to do multi. */
	int upnumvects; /* 0 if empty,  1 on bot level if not empty,
	 else # nonempty kids. */
	int *upnumeles; /* numeles[0] = # chgs on bot level, else
	 number of terms in kid's expansion. */
	double **upvects; /* vects[0] = chgs on bot level, else vectors
	 of kids' expansion terms. */
	int multisize; /* Number of terms in the expansion. */
	double *multi; /* Vector of multi coefficients. */
	double ***upmats; /* Matrices for chgs to multi or multi to multi.
	 upmats[i] is multisize x upnumeles[i]. */
	int *is_dummy; /* is_dummy[i] = TRUE => panel i is a dummy panel
	 used for elec field eval - omit from upward pass */
	int *is_dielec; /* is_dielec[i] = TRUE => panel i is on a surf
	 of type DIELEC or BOTH */

	/* Downward Pass variables. */
	int loc_exact; /* TRUE => do not build a local expansion */
	struct cube *lnext; /* Ptr to next cube on which to do local. */
	int downnumvects; /* Number of cubes in iteraction list. */
	int *downnumeles; /* # of eles in interact cube's expansion. */
	double **downvects; /* Vects of interact cube's expansion. */

	int localsize; /* Size of the local expansion */
	double *local; /* Vector of local field coefs */
	double ***downmats; /* Matrices for multi to chg, or multi to local
	 or local to local.  Downnumele x localsize. */

	struct cube **interList; /* explicit interaction list
	 - for fake dwnwd passes and eval pass */
	int interSize; /* number of elements in interList
	 - often != downnumvects nor evalnumvects */

	/* evaluation pass variables */
	struct cube *enext; /* pntr to next cube to evaluate */
	int evalnumvects; /* for exact = #in inter list, o.w. = 1 */
	int *evalnumeles; /* num of elements in inter list entry exp */
	double **evalvects; /* multi, local, or chgs of ilist entry */

	double *eval; /* vector of evaluation pnt voltages in cube */
	double ***evalmats; /* matrices for multi to potential, local to
	 potential or charge to potential */
	int *eval_isQ2P; /* is this matrix done directly?
	 For inductance stuff   12/92 */

	/* Direct portion variables. */
	struct cube *dnext; /* Ptr to next cube on which to do direct. */
	struct cube *pnext; /* Ptr to next cube on which to do precond. */
	struct cube *rpnext; /* Reverse ptr to next cube to do precond. */
	int dindex; /* Used to determine lower triang portion. */
	int directnumvects; /* Number of vects, self plus nbrs. */
	int *directnumeles; /* # of elements in the nbrs chg vect.
	 directnumeles[0] = numchgs in cube. */
	double **directq; /* Vecs of chg vecs, directq[0] this cube's. */
	double ***directmats; /* Potential Coeffs in cube and neighbors. */
	double ***precondmats; /* Precond Coeffs in cube and neighbors. */
	double **directlu; /* Decomposed cube potential Coefficients. */
	double **precond; /* Preconditioner. */
	double *prevectq; /* The charge vector for the preconditioner. */
	double *prevectp; /* The potential vector for preconditioner. */
	int presize; /* Size of the preconditioner. */
	int **nbr_is_dummy; /* Dummy vectors corresponding to directq's */
	//int **nbr_is_dummyo; /* BAP change */
	/* Cube structure variables. */
	charge **chgs; /* Array of charge ptrs. Only used lowest level. */
	struct cube **nbrs; /* Array of ptrs to nonemptry nearest neighbors. */
	int numnbrs; /* Number of nonempty neighbors. */
	struct cube **kids; /* Array of children ptrs. */
	int numkids; /* Number of kids. */
	struct cube *parent; /* Ptr to parent cube. */
	int parent_num;//BAP change
	int kid_num;//Bap change
	int path_open;//Bap change

};
typedef struct cube cube;

struct ssystem {
	int side; /* # cubes per side on lowest level. */
	int depth; /* # of levels of cubes. */
	int order; /* # of levels of cubes. */
	int num_cond; /* number of conductors */
	Name *cond_names; /* conductor name list */
	double perm_factor; /* overall scale factor for permittivities */
	double length; /* Length per cube on lowest level. */
	double minx, miny, minz; /* Coordinates of one corner of the domain. */
	int mul_maxq; /* max #panels in mul_exact cube */
	int mul_maxlq; /* max #panels in lowest level cube */
	int max_panel; /* max #panels in all cubes w/multipole */
	int up_size; /* total #panels in all cubes, incl dummies */
	int loc_maxq; /* max #evaluation points in loc_exact cube */
	int loc_maxlq; /* max #eval pnts in lowest level cube. */
	int max_eval_pnt; /* max #eval pnts in all cubes w/local exp */
	int eval_size; /* total #eval pnts in entire problem */
	double *q; /* The vector of lowest level charges. */
	double *p; /* The vector of lowest level potentials. */
	charge *panels; /* linked list of charge panels in problem */
	cube *****cubes; /* The array of cube pointers. */
	cube **multilist; /* Array of ptrs to first cube in linked list
	 of cubes to do multi at each level. */
	cube **locallist; /* Array of ptrs to first cube in linked list
	 of cubes to do local at each level. */
	cube *directlist; /* head of linked lst of low lev cubes w/chg */
	cube *precondlist; /* head of linked lst of precond blks. */
	cube *revprecondlist; /* reversed linked lst of precond blks. */
	int *is_dummy; /* is_dummy[i] = TRUE => panel i is a dummy */
	int *is_dielec; /* is_dielec[i] = TRUE => panel i on dielec */
	int DirectEval; /* TRUE - mulEval will do only Q2P products.
	 FALSE - mulEval will do all others. 12/92 */
	int real_depth; /* depth w/o fake levels for preconditioner*/
	int debug; /* print debug info */
	//start BAP
	int ***Parent_kids;
	int *number_level_groups;
	double ***Center;
	int **kids_index;
	double distance;
	charge **fill;
	int **kid_num;
	int ***child_index;
	int group_num;
	int **group_interaction_number;
	int observation;
	int unknown_num;
	//int *unknown_indices;
	int ***interactions;
	int parent_index;
	int **next_iteration;
	int *number_next_iteration;


	double ***U;
	double ***V;
	int count;
	int *k;
	int k_max;
	int el_max;
	double *direct_inductance;
	long long int memory;
	long long int memory_aa;
	//end Bap
/*
 int *****evalQ2Ps;
 int *****evalL2Ps;
 int *****evalM2Ps;
 int *****cntQ2Ps;
 int *****cntL2Ps;
 int *****cntM2Ps;
 */
};
typedef struct ssystem ssystem;

