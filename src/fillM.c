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


/* this is the new fillM    10/92 */

#include <string.h>
#include "induct.h"
#include "gp.h"

/* SRW */
void fillM(SYS*);
MELEMENT *make_mesh_from_path(SPATH*, int, SYS*);
int is_next_seg_in_gp(SPATH*, NODES*);
MELEMENT *insert_in_list(MELEMENT*, MELEMENT*);
NODES *getnode(int, seg_ptr);
void bad_seg_type(char*, seg_ptr);
MELEMENT *make_melement(int, FILAMENT*, int);
void add_to_external(PSEUDO_SEG*, int, int, SYS*);
int_list *make_int_list(int, int);
int_list *add_to_int_list(int_list*, int_list*);
int makeMlist(GROUNDPLANE*, MELEMENT**, Minfo*, int);
void fill_b(EXTERNAL*, CX*);
#ifdef SU0914
void extractYcol(CX **mat, CX *x0, EXTERNAL *extcol, EXTERNAL *ext_list,
		int write_currents);
#else
void extractYcol(CX**, CX*, EXTERNAL*, EXTERNAL*);
#endif // SU0914
char *get_a_name(PSEUDO_SEG*);
NODES *find_first_node(SPATH*);
void makegrids(SYS*, CX*, int, int);

/* this fills the kircoff's voltage law matrix (Mesh matrix) */
/* it maps a matrix of mesh currents to branch currents */
/* it might actually be what some think of as the transpose of M */
/* Here, M*Im = Ib  where Im are the mesh currents, and Ib the branch */
/* 6/92 I added Mlist which is a vector of linked lists to meshes. 
 This replaces M.  But I keep M around for checking things in matlab. */

/* much of what is commented out is obsolete stuff from an old idea
 for a preconditioner that never worked */
void fillM(SYS *indsys) {

	GROUNDPLANE *plane; /* CMS 6/2/92 */

	int mesh, k, minimeshes;
	MELEMENT **Mlist, *mend;
	Minfo *m_info;
	TREE *atree;
	PATHLIST *aplist;
	SPATH *apath;
	SEGMENT *seg;
	PSEUDO_SEG *pseg;

	Mlist = indsys->Mlist;
	m_info = indsys->m_info;
	mesh = 0;

	/* do all the loops due circuits in the graph */
	for (atree = indsys->trees; atree != NULL; atree = atree->next)
		for (aplist = atree->loops; aplist != NULL; aplist = aplist->next) {

			Mlist[mesh] = make_mesh_from_path(aplist->path, mesh, indsys);
			if (Mlist[mesh] != NULL)
				/* it's possible that coincident gp nodes cause no path at all */
				mesh++;

			/* unimplemented junk
			 make_unconstrained(&(m_info[mesh]), mesh);
			 mesh += make_many_meshes_from_path(aplist->path, Mlist, m_info, mesh,
			 indsys);
			 */
		}

	if (mesh > indsys->extra_meshes) {
		fprintf(stderr, "Internal Error: Bad estimate for extra_meshes\n");
		fprintf(stderr, "   One solution is to change FILS_PER_MESH to 1\n");
		exit(1);
	}

	minimeshes = mesh;

	/* this does all the meshes due to filaments within each segment */
	for (seg = indsys->segment; seg != NULL; seg = seg->next) {
		for (k = 1; k < seg->num_fils; k++, mesh++) {

			Mlist[mesh] = make_melement(seg->filaments[k - 1].filnumber,
					&seg->filaments[k - 1], 1);

			Mlist[mesh] = insert_in_list(
					make_melement(seg->filaments[k].filnumber,
							&seg->filaments[k], -1), Mlist[mesh]);
			/* unimplemented junk
			 make_unconstrained(&(m_info[mesh]),mesh);
			 */
		}
	}

	/* add all the lists to the groundplane */
	for (plane = indsys->planes; plane != NULL; plane = plane->next) {
		if (!is_nonuni_gp(plane))
			mesh += makeMlist(plane, &(Mlist[mesh]), &(m_info[mesh]), mesh);
		else
			mesh += make_nonuni_Mlist(plane, &(Mlist[mesh]));
	}

	/* For each tree mesh, pick one mini-mesh to be unconstrained and make it
	 unique (unimplemented)*/
	/*pick_unconstrained(Mlist, m_info, mesh, indsys->tree_meshes, minimeshes);*/

	if (mesh <= indsys->num_mesh)
		indsys->num_mesh = mesh;
	else {
		fprintf(stderr, "uh oh, mesh > num_mesh\n");
		exit(1);
	}
}

#ifdef SRW0814
/* SRW
 * A sort function for the MELEMENT list.  For long lists, it is much
 * faster to bulk sort than to keep the list ordered during insertion.
 */

static int mcmp(const void *a, const void *b) {
	MELEMENT *m1 = *(MELEMENT**) a;
	MELEMENT *m2 = *(MELEMENT**) b;
	return (m1->filindex - m2->filindex);
}

MELEMENT *msort(MELEMENT *m0) {
	int i, cnt = 0;
	MELEMENT **ary, *m;
	if (!m0 || !m0->mnext)
		return (m0);
	for (m = m0; m; m = m->mnext)
		cnt++;
	ary = (MELEMENT**) malloc(cnt * sizeof(MELEMENT*));
	cnt = 0;
	for (m = m0; m; m = m->mnext)
		ary[cnt++] = m;
	qsort(ary, cnt, sizeof(MELEMENT*), mcmp);
	cnt--;
	for (i = 0; i < cnt; i++)
		ary[i]->mnext = ary[i + 1];
	ary[cnt]->mnext = 0;
	m = ary[0];
	free(ary);
	return (m);
}

#endif

/* this takes a linked list of segments (path) and makes a row of the */
/* mesh matrix out of the filament[0]'s of each segment.  */
MELEMENT *make_mesh_from_path(SPATH *path, int mesh, SYS *indsys) {
	SPATH *selem, *temppath, *telem;
	MELEMENT *m1, *m2, *m3, *mlist;
	NODES *plusnode, *node, *plus2, *node0, *node1, *actualnode;
	int i, sign, sign2;
	SEGMENT *seg;
	PSEUDO_SEG *pseg, *pseg2;

	mlist = NULL;
	plusnode = find_first_node(path); /* which node starts the loop */
	for (selem = path, i = 0; selem != NULL; selem = selem->next, i++) {
		node0 = getnode(0, selem->seg); /* get original (not real) nodes */
		node1 = getnode(1, selem->seg);
		if (getrealnode(plusnode) == getrealnode(node0))
			sign = 1;
		else if (getrealnode(plusnode) == getrealnode(node1))
			sign = -1;
		else {
			fprintf(stderr,
					"make_mesh_from_path: segments don't connect at node %s\n",
					plusnode->name);
			exit(1);
		}

		if (selem->seg.type == NORMAL) {
			seg = (SEGMENT *) selem->seg.segp;

			m1 = make_melement(seg->filaments[0].filnumber, &seg->filaments[0],
					sign);
#ifdef SRW0814
			m1->mnext = mlist;
			mlist = m1;
#else
			mlist = insert_in_list(m1, mlist);
#endif
		} else if (selem->seg.type == PSEUDO) {
			pseg = (PSEUDO_SEG *) selem->seg.segp;

			if (pseg->type == EXTERNTYPE)
				add_to_external(pseg, mesh, sign, indsys);
			else if (pseg->type == GPTYPE) {
				while (is_next_seg_in_gp(selem, plusnode) == TRUE) {
					/* this is an ugly addition to make gp meshes smaller if two segs */
					/* come from the same gp.  It was commented out because of a bug  */
					/* which hopefully i've fixed, 3/96 */
					selem = selem->next;
					if (sign == 1) {
						plusnode = node1;
						node1 = getothernode(node1, selem->seg);
					} else {
						plusnode = node0;
						node0 = getothernode(node0, selem->seg);
					}
					if (indsys->opts->debug == ON)
						printf("Fixing extra long gp mesh in gp %s, mesh %d.\n",
								node0->gp->name, mesh);
				}
				if (sign == 1) {
					temppath = path_through_gp(node0, node1, node0->gp);
					plus2 = node0;
				} else {
					temppath = path_through_gp(node1, node0, node0->gp);
					plus2 = node1;
				}
				while (temppath != NULL) {
					telem = temppath;
					seg = (SEGMENT *) telem->seg.segp;
					if (is_nonuni_gp(node0->gp)) {
						/* must compare cell nodes, not actual seg nodes */
						if (plus2->gp_node == seg->node[0]->gp_node) {
							sign2 = 1;
							actualnode = seg->node[0];
						} else if (plus2->gp_node == seg->node[1]->gp_node) {
							sign2 = -1;
							actualnode = seg->node[1];
						} else {
							fprintf(stderr,
									"Hey, path_through_gp made nonconnected path!\n");
							exit(1);
						}
					} else {
						if (plus2 == seg->node[0])
							sign2 = 1;
						else if (plus2 == seg->node[1])
							sign2 = -1;
						else {
							fprintf(stderr,
									"Hey, path_through_gp made nonconnected path!\n");
							exit(1);
						}
						actualnode = plus2;
					}
					m1 = make_melement(seg->filaments[0].filnumber,
							&seg->filaments[0], sign2);
#ifdef SRW0814
					m1->mnext = mlist;
					mlist = m1;
#else
					mlist = insert_in_list(m1, mlist);
#endif
					plus2 = getothernode(actualnode, telem->seg);
					temppath = temppath->next;
					/*	  free(telem); */
				}
			} else {
				fprintf(stderr, "make_mesh_from_path: unknown pseudo_seg %d\n",
						pseg->type);
				exit(1);
			}
		} else {
			bad_seg_type("make_mesh_from_path", selem->seg);
		}
		plusnode = getothernode(plusnode, selem->seg);
	}

#ifdef SRW0814
	mlist = msort(mlist);
#endif

	if (mlist == NULL) {
		fprintf(stderr,
				"make_mesh_from_path: Possible loop of .external statements which is not allowed!\n");
		fprintf(stderr,
				" .external's (possibly equiv'ed nodes) which may make a loop:\n");
		for (selem = path, i = 0; selem != NULL; selem = selem->next, i++) {
			node0 = getnode(0, selem->seg); /* get original (not real) nodes */
			node1 = getnode(1, selem->seg);
			fprintf(stderr, "  %s  %s\n", node0->name, node1->name);
			/* the above could be made more useful by searching through
			 the master pseudo_nodes list for pseudo_nodes thatpoint to these */
		}
	}

	return mlist;
}

/* Check to see if the next segment is also from the same groundplane */
int is_next_seg_in_gp(SPATH *selem, NODES *plusnode) {
	PSEUDO_SEG *pseg2, *pseg1;
	NODES *othernode;

	if (selem->next != NULL && selem->next->seg.type == PSEUDO) {
		pseg1 = (PSEUDO_SEG *) selem->seg.segp;
		pseg2 = (PSEUDO_SEG *) selem->next->seg.segp;
		if (pseg2->type == GPTYPE && pseg1->node[0]->gp == pseg2->node[0]->gp) {
			othernode = getothernode(plusnode, selem->seg);
			if (othernode == pseg2->node[0] || othernode == pseg2->node[1])
				return TRUE;
			else if (plusnode == pseg2->node[0] || plusnode == pseg2->node[1])
				/* segs could be reversed in the path list.  added 3/96 */
				return TRUE;
		}
		/*
		 else
		 printf("Not an error: Two adjacent gp segs not in same gp %s, %s.\n",
		 pseg1->node[0]->gp->name, pseg2->node[0]->gp->name);
		 */
	}
	return FALSE;
}

/* this inserts melem into the linked list beginning with bigmhead. */
/* It inserts it to preserve increasing filindex order. */
MELEMENT *insert_in_list(MELEMENT *melem, MELEMENT *bigmhead) {
	MELEMENT *melem2;

	if (bigmhead == NULL)
		return melem;
	else {
		/* find where in the list to put melem */
		melem2 = bigmhead;
		if (melem2->filindex > melem->filindex) { /* put at beginning */
			melem->mnext = melem2;
			bigmhead = melem;
		} else { /* find its place in the middle of the list */
			while (melem2->mnext != NULL
					&& melem2->mnext->filindex < melem->filindex)
				melem2 = melem2->mnext;
			/* insert it in the middle */
			melem->mnext = melem2->mnext;
			melem2->mnext = melem;
		}

		return bigmhead;
	}
}

NODES *getnode(int number, seg_ptr seg) {
	if (seg.type == NORMAL)
		return ((SEGMENT *) seg.segp)->node[number];
	else if (seg.type == PSEUDO)
		return ((PSEUDO_SEG *) seg.segp)->node[number];
	else
		bad_seg_type("getnode", seg);
	return ((NODES*) NULL);
}

void bad_seg_type(char *name, seg_ptr seg) {
	fprintf(stderr, "%s: bad seg type: %d\n", name, seg.type);
	exit(1);
}

MELEMENT *make_melement(int filindex, FILAMENT *fil, int sign) {
	MELEMENT *melem;

	melem = (MELEMENT *) MattAlloc(1, sizeof(MELEMENT));
	melem->filindex = filindex;
	melem->fil = fil;
	melem->sign = sign;
	melem->mnext = NULL;

	return melem;
}

/* this keeps track of the meshes which contain the nodes of the */
/* .external statement.  This will have a voltage source in them */
/* and will need a 1 placed in the RHS corresponding to mesh number 'mesh' */
void add_to_external(PSEUDO_SEG *pseg, int mesh, int sign, SYS *indsys) {
	EXTERNAL *port;
	int realsign;

	port = indsys->externals;
	while (port != NULL && port->source != pseg)
		port = port->next;

	if (port == NULL) {
		fprintf(stderr, "Hey, supposed external segment isn't in list\n");
		exit(1);
	}

	realsign = -1 * sign; /* since this will be moved to RHS, change its sign */

	port->indices = add_to_int_list(make_int_list(mesh, realsign),
			port->indices);

}

int_list *make_int_list(int mesh, int sign) {
	int_list *elem;

	elem = (int_list *) Gmalloc(sizeof(int_list));
	elem->index = mesh;
	elem->sign = sign;
	elem->next = NULL;

	return elem;
}

int_list *add_to_int_list(int_list *int_elem, int_list *list) {
	int_elem->next = list;
	return int_elem;
}

/* makes the Mlist for the groundplane given a plane and parameters defining */
/* the current location of the overall Mlist.                                */
int makeMlist(GROUNDPLANE *plane, MELEMENT **pMlist, Minfo *pm_info, int mstart) {
	MELEMENT *melem;
	SEGMENT *seg = NULL;
	int counter, i, j, k;
	int signofelem;
	int a_hole;
	SEGMENT ***segs1 = plane->segs1;
	SEGMENT ***segs2 = plane->segs2;

	counter = 0;

	for (i = 0; i < plane->seg2; i++) {
		for (j = 0; j < plane->seg1; j++) {

			if (segs1[j][i] != NULL && segs2[j + 1][i] != NULL
					&& segs1[j][i + 1] != NULL && segs2[j][i] != NULL) {
				pMlist[counter] = NULL;

				for (k = 0; k < 4; k++) {
					switch (k) {
					case 0:
						seg = plane->segs1[j][i];
						signofelem = -1.0;
						break;
					case 1:
						seg = plane->segs2[j + 1][i];
						signofelem = -1.0;
						break;
					case 2:
						seg = plane->segs1[j][i + 1];
						signofelem = 1.0;
						break;
					case 3:
						seg = plane->segs2[j][i];
						signofelem = 1.0;
						break;
					}

					melem = make_melement(seg->filaments[0].filnumber,
							&seg->filaments[0], signofelem);
					pMlist[counter] = insert_in_list(melem, pMlist[counter]);

				}
				/* unimplemented junk
				 make_unconstrained(&(pm_info[counter]), mstart+counter);
				 */
				counter++;
			}
		}
	}

	if (counter > plane->numesh) {
		printf("something wrong with meshes, numesh != counter \n");
		exit(1);
	}

	return counter;

}

void fill_b(EXTERNAL *ext, CX *b) {
	int_list *elem;

	for (elem = ext->indices; elem != NULL; elem = elem->next)
		b[elem->index].real = elem->sign;
}

#ifdef SU0914
void extractYcol(CX **mat, CX *x0, EXTERNAL *extcol, EXTERNAL *ext_list,
		int write_currents) {
	EXTERNAL *ext;
	int_list *elem;
	CX sum, tmp;

	if (write_currents)
		start_write_currents();

	for (ext = ext_list; ext != NULL; ext = ext->next) {
		sum = CXZERO;
		/* for each mesh that contains this voltage source */
		for (elem = ext->indices; elem != NULL; elem = elem->next) {
			cx_scalar_mult(tmp, elem->sign, x0[elem->index]);
			cx_add(sum, sum, tmp);
		}
		if (write_currents)
			place_current_output(extcol, ext, sum);
		printf("Port %s: %g%s%gj\n", ext->portname, sum.real,
				((sum.imag >= 0) ? "+" : ""), sum.imag);
		mat[ext->Yindex][extcol->col_Yindex] = sum;
		
	}

	


	if (write_currents)
		end_write_currents();
}
#else
void extractYcol(CX **mat, CX *x0, EXTERNAL *extcol, EXTERNAL *ext_list) {
	EXTERNAL *ext;
	int_list *elem;
	CX sum, tmp;
	for (ext = ext_list; ext != NULL; ext = ext->next) {
		sum = CXZERO;
		/* for each mesh that contains this voltage source */
		for (elem = ext->indices; elem != NULL; elem = elem->next) {
			cx_scalar_mult(tmp, elem->sign, x0[elem->index]);
			cx_add(sum, sum, tmp);
		}
		mat[ext->Yindex][extcol->col_Yindex] = sum;
	}
}
#endif // SU0914

char *get_a_name(PSEUDO_SEG *pseg) {
	return pseg->node[0]->name;
}

/* we wish to find the first node in a path which is the node of 
 the first segment which is not connected to the second segment */
NODES *find_first_node(SPATH *path) {
	NODES *node0, *node1;
	int node1_in_middle;
	int node0_in_middle;

	node0 = getnode(0, path->seg);
	node1 = getnode(1, path->seg);
	if (path->next == NULL) /* there is no other segment */
		return node0;

#if 1==0
	the old way

	replaced for the sake of fixing extra long gp meshes, we must
	handle a two segment loop more carefully /*

	 if (getrealnode(node0) == getrealnode(getnode(0, path->next->seg))
	 || getrealnode(node0) == getrealnode(getnode(1, path->next->seg)) )
	 /* node0 is connected to the next segment, so start with node 1 */
	return node1;
	else if (getrealnode(node1) == getrealnode(getnode(0, path->next->seg))
			|| getrealnode(node1) == getrealnode(getnode(1, path->next->seg)) )
	/* node0 is connected to the next segment, so start with node 1 */
	return node0;
	else {
		fprintf(stderr, "find_first_node: first seg not connected to second\n");
		exit(1);
	}
#endif

	/* is node0 connected to the next segment? */
	node0_in_middle = (getrealnode(node0)
			== getrealnode(getnode(0, path->next->seg))
			|| getrealnode(node0) == getrealnode(getnode(1, path->next->seg)));

	node1_in_middle = (getrealnode(node1)
			== getrealnode(getnode(0, path->next->seg))
			|| getrealnode(node1) == getrealnode(getnode(1, path->next->seg)));

	/* return the node that isn't connecting the first and second segs */
	if (node1_in_middle && !node0_in_middle)
		return node0;
	else if (node0_in_middle && !node1_in_middle)
		return node1;
	else if (node0_in_middle && node1_in_middle) {
		/* this is a two segment loop, so it doesn't matter which
		 we return.  But if these are two groundplane segments, perhaps
		 this needs to be shortened to one segment and in order to do so
		 we must determine the connectivity based on original, not real, node
		 name */
		if (node0 == getnode(0, path->next->seg)
				|| node0 == getnode(1, path->next->seg))
			return node1;
		else if (node1 == getnode(0, path->next->seg)
				|| node1 == getnode(1, path->next->seg))
			return node0;
		else
			/* it doesn't matter, both are equiv'd */
			return node0;
	} else {
		fprintf(stderr, "find_first_node: first seg not connected to second\n");
		exit(1);
	}

}

void makegrids(SYS *indsys, CX *Im, int column, int freq_num) {
	static CX *Ib = NULL, current;
	int fils, meshes;
	static CX **out1 = NULL;
	static CX **out2 = NULL;
	static int maxdir1 = 0, maxdir2 = 0;
	int dir1, dir2, num, i, j;
	MELEMENT *mtemp;
	GROUNDPLANE *p;
	FILE *fp, *fpreal, *fpimag, *fpmag;
	static char *fname, *tempstr;
	SEGMENT *seg;
	FILAMENT *fil;
	double xv, yv, zv, x, y, z, magcur;

	fils = indsys->num_fils;
	meshes = indsys->num_mesh;

	if (Ib == NULL) {
		Ib = (CX *) MattAlloc(fils, sizeof(CX));
		fname = malloc(100 * sizeof(char));
		tempstr = malloc(100 * sizeof(char));
	}

	/* do  Ib = Mtrans*Im */
	for (i = 0; i < fils; i++) {
		Ib[i] = CXZERO;
		for (mtemp = indsys->Mtrans[i]; mtemp != NULL; mtemp = mtemp->mnext) {
			if (mtemp->sign == 1)
				cx_add(Ib[i], Ib[i], Im[mtemp->filindex]);
			else
				cx_sub(Ib[i], Ib[i], Im[mtemp->filindex]);
		}
	}

	/* SRW -- stole Enrico's mod from FastHenry2 */
	// Enrico, added column / frequency information
	printf("saving to Jreal%s%d_%d.mat, Jimag%s%d_%d.mat, Jmag%s%d_%d.mat\n",
			indsys->opts->suffix, column + 1, freq_num, indsys->opts->suffix,
			column + 1, freq_num, indsys->opts->suffix, column + 1, freq_num);

	sprintf(fname, "Jreal%s%d_%d.mat", indsys->opts->suffix, column + 1,
			freq_num);

	/*
	 printf("saving to Jreal%s.mat, Jimag%s.mat, Jmag%s.mat\n",
	 indsys->opts->suffix, 
	 indsys->opts->suffix, 
	 indsys->opts->suffix);

	 sprintf(fname, "Jreal%s.mat",indsys->opts->suffix);
	 */

	/* SRW -- this is ascii data */
	fpreal = fopen(fname, "w");
	if (fpreal == NULL) {
		printf("couldn't open file %s\n", fname);
		exit(1);
	}
	/*  fprintf(fpreal, "$ DATA=VECTOR\n");*/

	sprintf(fname, "Jimag%s%d_%d.mat", indsys->opts->suffix, column + 1,
			freq_num);
	/*
	 sprintf(fname, "Jimag%s.mat",indsys->opts->suffix);
	 */
	/* SRW -- this is ascii data */
	fpimag = fopen(fname, "w");
	if (fpimag == NULL) {
		printf("couldn't open file %s\n", fname);
		exit(1);
	}
	/*  fprintf(fpimag, "$ DATA=VECTOR\n");*/

	sprintf(fname, "Jmag%s%d_%d.mat", indsys->opts->suffix, column + 1,
			freq_num);
	/*
	 sprintf(fname, "Jmag%s.mat",indsys->opts->suffix);
	 */
	/* SRW -- this is ascii data */
	fpmag = fopen(fname, "w");
	if (fpmag == NULL) {
		printf("couldn't open file %s\n", fname);
		exit(1);
	}
	/*  fprintf(fpmag, "$ DATA=VECTOR\n");*/

	for (seg = indsys->segment; seg != NULL; seg = seg->next)
		for (i = 0; i < seg->num_fils; i++) {
			fil = &(seg->filaments[i]);
			current = Ib[fil->filnumber];
			magcur = cx_abs(current);
			xv = fil->lenvect[XX] / fil->length / fil->area;
			yv = fil->lenvect[YY] / fil->length / fil->area;
			zv = fil->lenvect[ZZ] / fil->length / fil->area;
			x = fil->x[0];
			y = fil->y[0];
			z = fil->z[0];
			fprintf(fpreal, "%lg %lg %lg  %lg %lg %lg\n", x, y, z,
					xv * current.real, yv * current.real, zv * current.real);
			fprintf(fpimag, "%lg %lg %lg  %lg %lg %lg\n", x, y, z,
					xv * current.imag, yv * current.imag, zv * current.imag);
			fprintf(fpmag, "%lg %lg %lg  %lg %lg %lg\n", x, y, z, xv * magcur,
					yv * magcur, zv * magcur);
		}

	fclose(fpreal);
	fclose(fpimag);
	fclose(fpmag);

	if (indsys->num_planes == 0)
		return;

	printf("saving to file Grid%s%d_%d...\n", indsys->opts->suffix, column + 1,
			freq_num);
	sprintf(fname, "Grid%s%d_%d.mat", indsys->opts->suffix, column + 1,
			freq_num);

	/* SRW -- this is binary data */
	fp = fopen(fname, "wb");
	if (fp == NULL) {
		printf("couldn't open file %s\n", fname);
		exit(1);
	}

	for (num = 0, p = indsys->planes; p != NULL; p = p->next, num++) {
		if (is_nonuni_gp(p))
			dump_nonuni_plane_currents(p->nonuni, Ib, fp);
		else {
			dir1 = p->seg1 + 1;
			dir2 = p->seg2 + 1;

			if (dir1 > maxdir1 || dir2 > maxdir2) {
				out1 = (CX **) MatrixAlloc(dir2 + 10, dir1 + 10, sizeof(CX));
				out2 = (CX **) MatrixAlloc(dir2 + 10, dir1 + 10, sizeof(CX));
				maxdir1 = dir1 + 10;
				maxdir2 = dir2 + 10;
			}

			for (i = 0; i < dir2; i++)
				for (j = 0; j < dir1; j++) {
					/* do direction 1 */
					if (j != dir1 - 1 && p->segs1[j][i] != NULL) {
						out1[i][j] = Ib[p->segs1[j][i]->filaments[0].filnumber];
						if (p->segs1[j][i]->node[0] != p->pnodes[j][i]) {
							printf("You goofed 1\n");
						}
					} else
						out1[i][j] = CXZERO;

					/* do direction 2 */
					if (i != dir2 - 1 && p->segs2[j][i] != NULL) {
						out2[i][j] = Ib[p->segs2[j][i]->filaments[0].filnumber];
						if (p->segs2[j][i]->node[0] != p->pnodes[j][i]) {
							printf("You goofed 2\n");
						}
					} else
						out2[i][j] = CXZERO;
				}

			printf("saving grid1%s...\n", p->name);
			strcpy(fname, "grid1");
			sprintf(tempstr, "%s", p->name);
			strcat(fname, tempstr);

			savecmplx(fp, fname, out1, dir2, dir1);

			printf("saving grid2%s...\n", p->name);
			strcpy(fname, "grid2");
			sprintf(tempstr, "%s", p->name);
			strcat(fname, tempstr);

			savecmplx(fp, fname, out2, dir2, dir1);
		}
	}
	fclose(fp);
}
/*------------------------------------------------------------------------*/

#if 1==0

The following is code that was not implemented and never used and is now
obsolete.

/* this takes a linked list of segments (path) and makes many rows of the */
/* mesh matrix out of the filament[0]'s of each segment.  */
void make_many_meshes_from_path(SPATH *path, MELEMENT **Mlist, Minfo *m_info,
		int mstart, SYS *indsys)
/* int mstart;  the mesh index at which to start creating meshes. */
{
	SPATH *selem, *temppath, *telem, *minipath;
	MELEMENT *m1, *m2, *m3, *mlist;
	NODES *plusnode, *node, *plus2, *node0, *node1;
	int i, sign, sign2;
	SEGMENT *seg;
	PSEUDO_SEG *pseg, *pseg2;
	int mesh, seg_count;
	int physically_connected;
	PSEUDO_SEG *extern_seg;
	int extern_count, extern_sign;

	mlist = NULL;
	mesh = 0;
	minipath = NULL;
	seg_count = 0;
	extern_count = 0;

	plusnode = find_first_node(path); /* which node starts the loop */
	for(selem = path, i = 0; selem != NULL; selem = selem->next, i++) {
		node0 = getnode(0, selem->seg);
		node1 = getnode(1, selem->seg);
		if (getrealnode(plusnode) == getrealnode(node0))
		sign = 1;
		else if (getrealnode(plusnode) == getrealnode(node1))
		sign = -1;
		else {
			fprintf(stderr, "make_mesh_from_path: segments don't connect at node %s\n", plusnode->name);
			exit(1);
		}

		if (selem->seg.type == NORMAL) {
			physically_connected = seg_count != 0 && (sign==1 && plusnode == node0
					|| sign==-1 && plusnode==node1);
			if (physically_connected) {
				minipath = add_seg_to_list(selem->seg, minipath);
				seg_count++;
			}
			if (!physically_connected || seg_count == FILS_PER_MESH) {
				Mlist[mstart + mesh] = make_mesh_from_path(minipath, mstart+mesh,
						indsys);
				mesh++;
				reset_vars(&seg_count, &minipath);
			}
		}
		else if (selem->seg.type == PSEUDO) {
			pseg = (PSEUDO_SEG *)selem->seg.segp;

			if (pseg->type == EXTERNTYPE) {
				extern_count++;
				extern_seg = pseg;
				extern_sign = sign;
				/*	add_to_external(pseg, mesh, sign, indsys); */
			}
			else if (pseg->type == GPTYPE) {
				/* flush any mesh which isn't finished yet */
				if (seg_count != 0) {
					Mlist[mstart + mesh] = make_mesh_from_path(minipath, mstart+mesh,
							indsys);
					mesh++;
					reset_vars(&seg_count, &minipath);
				}
				while( is_next_seg_in_gp(selem) == TRUE && 1==0) {
					/* this is an ugly addition to make gp meshes smaller if two segs */
					/* come from the same gp */
					selem = selem->next;
					if (sign == 1) {
						plusnode = node1;
						node1 = getothernode(node1, selem->seg);
					}
					else {
						plusnode = node0;
						node0 = getothernode(node0, selem->seg);
					}
					if (indsys->opts->debug == ON)
					printf("Fixing extra long gp mesh in gp %s, mesh %d.\n",
							node0->gp->name, mesh);
				}
				if (sign == 1) {
					temppath = path_through_gp(node0, node1, node0->gp);
					plus2 = node0;
				}
				else {
					temppath = path_through_gp(node1, node0, node0->gp);
					plus2 = node1;
				}
				while(temppath != NULL) {
					telem = temppath;
					seg = (SEGMENT *)telem->seg.segp;
					if (plus2 == seg->node[0])
					sign2 = 1;
					else if (plus2 == seg->node[1])
					sign2 = -1;
					else {
						fprintf(stderr, "Hey, path_through_gp made nonconnected path!\n");
						exit(1);
					}
					minipath = add_seg_to_list(selem->seg, minipath);
					seg_count++;
					if (seg_count == FILS_PER_MESH) {
						Mlist[mstart + mesh] = make_mesh_from_path(minipath, mstart+mesh,
								indsys);
						mesh++;
						reset_vars(&seg_count, &minipath);
					}

					plus2 = getothernode(plus2, telem->seg);
					temppath = temppath->next;
				}
			}
			else {
				fprintf(stderr, "make_mesh_from_path: unknown pseudo_seg %d\n",
						pseg->type);
				exit(1);
			}
		}
		else {
			bad_seg_type("make_mesh_from_path", selem->seg);
		}
		plusnode = getothernode(plusnode, selem->seg);
	}

	if (seg_count == FILS_PER_MESH) {
		Mlist[mstart + mesh] = make_mesh_from_path(minipath, mstart+mesh,
				indsys);
		mesh++;
		reset_vars(&seg_count, &minipath);
	}

	if (extern_count > 0) {
		if (extern_count == 1)
		add_to_external(extern_seg, mstart, extern_sign, indsys);
		else {
			fprintf(stderr,"Internal Err: More than one voltage source in a mesh\n");
			exit(1);
		}
	}

#if 1==0
	/* only the first mini mesh is unconstrained */
	make_m_info(&(m_info[0]), UNCONSTRAINED, mstart, mstart, mstart, mesh);

	/* make the second one constrained to the first */
	if (mesh > 1) {
		make_m_info(&(m_info[1]), CONSTRAINED, mstart+1, mstart, mstart, mesh);
	}

	/* The first constrained mesh above will have a current for these */
	/* i.e., they are constrained to mstart+1 */
#endif

	for(i = 0; i < mesh; i++) {
		make_m_info(&(m_info[i]), CONSTRAINED, mstart+i, -1, mstart,mesh, 0);
	}

	if (mesh == 0) {
		fprintf(stderr, "make_many_meshes: No meshes made!\n");
		exit(1);
	}

	return mesh;
}

void reset_vars(nt *seg_count, SPATH **minipath)
{
	*seg_count = 0;
	free_spath(*minipath);
	*minipath = NULL;
}

void make_unconstrained(Minfo *pm_info, int mesh)
{
	make_m_info(pm_info, UNCONSTRAINED, mesh, -1, -1, -1, -1);
}

void make_m_info(Minfo *pm_info, int type, int mesh_num, int constr_mesh,
		int first, int num_meshes, int other_mesh)
{
	pm_info->type = type;
	pm_info->mesh_num = mesh_num;
	pm_info->constraining_mesh = constr_mesh;
	pm_info->first = first;
	pm_info->num_meshes = num_meshes;
	pm_info->other_mesh = other_mesh;
}

void pick_unconstrained(MELEMENT **Mlist, Minfo *m_info, int total_meshes,
		int big_meshes, int num_meshes)
{
	Minfo **m_undone, *m_one;
	int i, j, quit;
	int counter;
	int last_undone, undone;
	int first /*, num_meshes*/;

	counter = 0;
	m_undone = (Minfo **)Gmalloc(big_meshes*sizeof(Minfo *));

	/* This goes through all the big meshes and picks one representative mesh*/
	for(i = 0; i < num_mesh; i++) {
		m_undone[counter++] = &(m_info[i]);
		i += m_info[i].num_meshes;
	}

	if (counter != big_meshes) {
		fprintf(stderr, "Error getting representative meshes\n");
		exit(1);
	}

	undone = big_meshes;
	last_undone = undone + 1;

	while(undone < last_undone) {
		last_undone = undone;

		/* counts duplicates and saves them in m_undone[i]->other_mesh */
		count_duplicates(m_undone, undone, Mlist, m_info);

		for(i = 0; i < undone; i++) {
			first = m_undone[i]->first;
			num_meshes = m_undone[i]->num_meshes;
			quit = FALSE;
			for(j = first; j < first + num_meshes && quit == FALSE; j++) {
				if (!is_duplicated(m_info[j])) {
					/* This mesh is unique among minimeshes*/
					if (is_globally_unique(m_info[j], num_mesh, Mlist, total_meshes)) {
						quit = TRUE;

						choose_this_mesh(m_undone[i], j, m_info);
					}
				}
			}
		}

		undone = get_undone(m_undone, undone);
	}

	/* sort in order of increasing number of minimeshes */
	/* This is a rock sort.  Biggest element falls to the bottom */
	for(i = 0; i < undone - 1; i++)
	for(j = 1; j < undone - i; j++)
	if (m_undone[j-1]->num_meshes > m_undone[j]->num_meshes) {
		m_one = m_undone[j-1];
		m_undone[j-1] = m_undone[j];
		m_undone[j] = m_one;
	}

	/* Go through all the undone meshes and pick the first available mesh */
	for(i = 0; i < undone; i++) {
		first = m_undone[i]->first;
		num_meshes = m_undone[i]->num_meshes;
		quit = FALSE;
		for(j = first; j < first + num_meshes && quit == FALSE; j++) {
			if (is_locally_unique(m_undone, j, Mlist)) {
				/* This mesh is unique among other selected */
				if (is_globally_unique(m_info[j], &(Mlist[num_mesh]),
								total_meshes - num_mesh)) {
					quit = TRUE;

					choose_this_mesh(m_undone[i], j, m_info);
				}
			}
		}
		if (quit == FALSE) {
			fprintf(stderr, "screwed up!\n");
			exit(1);
		}
	}

}

void count_duplicates(Minfo **m_undone, int undone, MELEMENT **Mlist,
		Minfo *m_info)
{
	int i,j,k,m;
	Minfo *m_one;
	int first_s, num_meshes_s, first_t, num_meshes_t;

	for(i = 0; i < undone; i++) {
		first_s = m_undone[i]->first;
		num_meshes_s = m_undone[i]->num_meshes;
		for(j = first_s; j < first_s + num_meshes_s; j++) {
			for(k = 0; k < undone; k++) {
				first_t = m_undone[i]->first;
				num_meshes_t = m_undone[i]->num_meshes;
				for(m = first_t; m < first_t + num_meshes_t; m++)
				if(mesh_comp(Mlist[j], Mlist[m]) == TRUE) {
					m_info[j].other_mesh += 1;
					if (j != m)
					m_info[m].other_mesh += 1;
				}
			}
		}
		if (m_info[j].other_mesh == 0) {
			fprintf(stderr, "Huh?  mesh isn't identical to itself?\n");
			exit(1);
		}
	}
}

void mesh_comp(MELEMENT *m1, MELEMENT *m2)
{
	int same;

	same = TRUE;
	while(m1 != NULL && m2 != NULL && same == TRUE) {
		if (m1->filindex != m2->filindex)
		same = FALSE;
		else {
			m1 = m1->next;
			m2 = m2->next;
		}

		if (m1 != NULL || m2 != NULL)
		same = FALSE;

		return same;
	}

	int is_duplicated(Minfo m_one)
	{
		return (m_one.other_mesh != 1);
	}

	int is_globally_unique(Minfo m_one, int minimeshes, MELEMENT **Mlist,
			int num_mesh)
	{
		int unique = TRUE;
		int i,j;

		for(i = minimeshes; i < num_mesh && unique == TRUE; i++)
		if (mesh_comp(Mlist[m_one.mesh_num], Mlist[i]) == TRUE)
		unique = FALSE;

		return unique;
	}

	void choose_this_mesh(Minfo *m_begin, int constraining_mesh, Minfo *m_info)
	{
		int other_mesh = -1;
		int num_dups = 1000000;

		m_info[constraining_mesh] = UNCONSTRAINED;

		/* find a good other_mesh */
		for(i = m_begin->first; i < m_begin->first + m_begin->num_meshes; i++)
		if (i != constraining_mesh)
		if (m_info[i].other_mesh < num_dups) {
			other_mesh = i;
			num_dups = m_info[i].other_mesh;
		}

		if (m_begin->num_meshes != 1) {
			if (other_mesh == -1) {
				fprintf(stderr, "Huh?  other_mesh == -1?\n");
				exit(1);
			}
			for(i = m_begin->first; i < m_begin->first + m_begin->num_meshes; i++) {
				m_info[i].constraining_mesh = constraining_mesh;
				m_info[i].other_mesh = other_mesh;
			}
		}
	}

	/* get all the remaining undone big meshes based on constraing_mesh still -1*/
	int get_undone(Minfo **m_undone, int undone)
	{
		int i, new_undone;

		new_undone = 0;
		for(i = 0; i < undone; i++)
		if (m_undone[i]->constraining_mesh == -1)
		m_undone[new_undone++] = m_undone[i];

		return new_undone;
	}

	int is_locally_unique(Minfo **m_undone, int mini_mesh, MELEMENT **Mlist)
	{
		int i, k;
		int unique = TRUE;

		for(i = 0; i < j && unique == TRUE; i++)
		if (mesh_comp(Mlist[j], Mlist[m_undone[i]->constraining_mesh]) == TRUE)
		unique = FALSE;

		return unique;
	}

#endif

