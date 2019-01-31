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




/* this is 'fastcap.c' */
#include "induct.h"

/*Here mulMatDirect is still used as the preconditioner requires this.
For further optimization this can be replaced using the same integration method used by the MLACA - SVD.
All three mulMatUp, mulMatDown and mulMatEval are only required when the MLFMA approach is used.*/

/* SRW */
ssystem *SetupMulti(charge*, SYS*);

double aanalitic(FILAMENT *, FILAMENT *, int, int, int);//bap change

ssystem *SetupMulti(charge *chglist, SYS *indsys) {
	int ttliter, i, j, num_cond;
	charge *nq;
	ssystem *sys;
	double **capmat, dirtimesav, mulsetup, initalltime, ttlsetup, ttlsolve;
	double relperm;
	int autmom, autlev, numMom, numLev = 0;

	int pseudo_lev; /* level for breaking up only */

	extern int fulldirops, fullPqops;
	extern int num_dummy_panels, num_dielec_panels;
	extern int num_both_panels, num_cond_panels, up_size, eval_size;
	extern char *title, *ps_file_base;
	extern long memcount;
	extern double prectime, conjtime, dirtime, multime, uptime, downtime;
	extern double evaltime, lutime, fullsoltime, prsetime;
	extern char *kill_name_list;

	Name *name_list;

#if DUMPPS == ON || DUMPPS == ALL
	char filename[BUFSIZ];
#endif

#if CAPVEW == ON
	extern char **argvals;
	extern int argcnt, m_, q_, dd_;
	extern double ***axes;
#endif

	/* initialize memory and time counters, etc. */
	fulldirops = fullPqops = 0;
	prectime = conjtime = dirtime = multime = uptime = downtime = 0.0;
	evaltime = lutime = fullsoltime = mulsetup = 0.0;
	/*  memcount = 0;   Do this in induct.c */
	/*  CALLOC(title, BUFSIZ, char, ON, AMSC); */

	/* initialize defaults, etc */
	autmom = autlev = ON;
	relperm = 1.0;

#if CAPVEW == ON
	/*  argvals = argv;
	argcnt = argc;
	*/
	CALLOC(axes, 10, double **, ON, AMSC);
	for (i = 0; i < 10; i++) {
		CALLOC(axes[i], 2, double *, ON, AMSC);
		for (j = 0; j < 2; j++) {
			CALLOC(axes[i][j], 3, double, ON, AMSC);
		}
	}
#endif

	/* get the list of all panels in the problem */
	/* - many command line parameters having to do with the postscript
	file dumping interface are passed back via globals (see mulGlobal.c) */
	/*  chglist = input_problem(argv, argc, &autmom, &autlev, &relperm,
	&numMom, &numLev, &name_list, &num_cond);
	*/
	if (autmom == ON)
		numMom = DEFORD;
	else
		numMom = DEFORD;

#if DIRSOL == ON || EXPGCR == ON
	/*fprintf(stderr, "DIRSOL and EXPGCR compile options not implemented\n");
	exit(1);*/
	*numLev = 0; /* put all the charges in first cube */
	*autlev = OFF;
#endif

#if CAPVEW == ON
	/* if no fastcap run is to be done, just dump the psfile */
	/*
	if(m_) {
	if(!q_) get_ps_file_base(argv, argc);
	dump_ps_geometry(chglist, NULL, 0, dd_);
	exit(0);
	}
	*/
#endif

	/* the following is added for FastHenry to alter settings for mulInit()*/
	//indsys->opts->level = 7;//BAP change
	autlev = OFF;
	numMom = indsys->opts->order;
	if (indsys->opts->level == AUTO)
		autlev = ON;
	else
		numLev = indsys->opts->level;

	/*
	if (numLev == 0) {
	printf("Level to break up segments:");
	scanf("%d",&pseudo_lev);
	}
	else
	*/
	pseudo_lev = numLev;

	/*  starttimer; */
	sys = mulInit(autlev, numLev, numMom, chglist, indsys, pseudo_lev);
	/*Set up cubes, charges.*/
	/* stoptimer */
	initalltime = dtime;
	numLev = sys->depth;

	/*
	sys->num_cond = num_cond;
	sys->cond_names = name_list;
	*/

	//fprintf(stdout, "\nMultipole Summary\n");

#if CMDDAT == ON
	
	//fprintf(stdout, "  Expansion order: %d\n", numMom);
	if (sys->depth == sys->real_depth)
	fprintf(stdout, "  Number of partitioning levels: %d\n", numLev);
	else {
	fprintf(stdout, "  Number of partitioning levels: %d\n",
	sys->real_depth);
	fprintf(stdout, "  Number of levels for Preconditioner: %d\n", numLev);
	}
	
	/*  fprintf(stdout, "  Overall permittivity factor: %.3g\n", relperm); */
#endif

	/* Figure out number of panels and conductors. */
	eval_size = up_size = num_dummy_panels = num_dielec_panels = 0;
	num_both_panels = num_cond_panels = 0;
	for (nq = chglist; nq != NULL; nq = nq->next) {
		if (nq->dummy)
			num_dummy_panels++;
		else if (nq->surf->type == CONDTR) {
			num_cond_panels++;
		}
		else if (nq->surf->type == DIELEC)
			num_dielec_panels++;
		else if (nq->surf->type == BOTH)
			num_both_panels++;
	}
	up_size = num_cond_panels + num_both_panels + num_dielec_panels;
	eval_size = up_size + num_dummy_panels;

#if DISSRF == OFF
	/*  fprintf(stdout, "Title: `%s'\n", title); */
#endif
	fprintf(stdout, "  Total number of filaments: %d\n", up_size);
#if 1==0
	fprintf(stdout, "    Number of conductor panels: %d\n", num_cond_panels);
	fprintf(stdout, "    Number of dielectric interface panels: %d\n",
		num_dielec_panels);
	fprintf(stdout,
		"    Number of thin conductor on dielectric interface panels: %d\n",
		num_both_panels);
	/*fprintf(stdout, "  Number of extra evaluation points: %d\n",
	num_dummy_panels);*/
	fprintf(stdout, "  Number of conductors: %d\n", num_cond);
#endif

	/*********************************************************************/
	/* exit now if we only want to do multipole setup and never evaluate */
	if (indsys->opts->mat_vect_prod == DIRECT) {
		/* fprintf(stdout, "  Inductance mem allocated: %7.d kilobytes\n",
		memIND/1024); */
		//Start Bap
		//ssystem *syss = sys;

		//int track = 0;
		//mulMatDirect(syss);//Bap change
		/*
		sys->number_groups_level = (int *)malloc(sys->depth * sizeof(int));
		sys->number_of_interactions_groups = (int **)malloc(sys->depth * sizeof(int*));
		sys->number_els_sys = (int **)malloc(sys->depth * sizeof(int*));
		sys->number_els_sys = (int ***)malloc(sys->depth * sizeof(int**));
		sys->k = (int *)malloc(syss->count * sizeof(int));
		sys->U = (int ***)malloc(syss->count * sizeof(int**));
		sys->V = (int ***)malloc(syss->count * sizeof(int**));
		for (int i = 0; i < sys->depth -1; i++)
		{
		sys->number_of_interactions_groups[i] = (int *)malloc(syss->number_groups_level[i] * sizeof(int));
		sys->number_els_sys[i] = (int *)malloc(syss->number_groups_level[i] * sizeof(int));
		sys->number_groups_level[i] = syss->number_groups_level[i];
		for (int j = 0; j < syss->number_groups_level[i]; j++)
		{
		sys->number_of_interactions_groups[i][j] = sys->number_of_interactions_groups[i][j];
		sys->number_els_sys[i][j] = syss->number_els_sys[i][j];
		for (int y = 0; y < syss->number_of_interactions_groups[i][j]; y++)
		{

		sys->U[track] = (int **)malloc(syss->number_els_sys[i + 1][j] * sizeof(int*));
		sys->V[track] = (int **)malloc(syss->k[track] * sizeof(int*));


		for (int o = 0; o < syss->number_els_sys[i+1][j]; o++)
		{
		sys->U[track][o] = (int *)malloc(syss->k[track] * sizeof(int));
		for (int kk = 0; kk < syss->k[track]; kk++)
		{
		sys->U[track][o][kk] = syss->U[track][o][kk];
		}
		}
		for (int kk = 0; kk < syss->k[track]; kk++)
		{
		sys->V[track][kk] = (int *)malloc(syss->number_els_sys[i + 1][syss->interact[i][j][y]] * sizeof(int));
		for (int o = 0; o < syss->number_els_sys[i+1][syss->interact[i][j][y]]; o++)
		{

		sys->V[track][kk][o] = syss->V[track][kk][o];

		}
		}
		sys->k[track] = syss->k[track];
		++track;
		//memcpy(sys->U, syss->U, ???);
		//memcpy(sys->V, syss->V, ???);
		}
		}
		}
		*/

		/*
		for (int i = 0; i < indsys->num_fils; i++)
		{
		for (int j = 0; j < indsys->num_fils; j++)
		{
		//	sys->matrix[i][j] = syss->matrix[i][j];
		}
		}
		*/
		//End Bap
		return sys;
	}
	/*********************************************************************/

#if NAMDAT == ON
	dumpCondNames(stdout, name_list);
#endif

	if (num_both_panels > 0) {
		fprintf(stderr,
			"Thin cond panels on dielectric interfaces not supported\n");
		exit(1);
	}

#if CKCLST == ON
	fprintf(stdout, "Checking panels...");
	if (has_duplicate_panels(stdout, chglist)) {
		fprintf(stdout, "charge list has duplicates\n");
		exit(1);
	}
	fprintf(stdout, "no duplicates\n");
#endif

#if MULDAT == ON
	dumpMulSet(sys, numLev, numMom);
#endif
	fflush(stdout);

	/* starttimer; */
	mulMultiAlloc(MAX(sys->max_eval_pnt, sys->max_panel), numMom, sys->depth);
	/* stoptimer */
	initalltime += dtime; /* save initial allocation time */

#if DUMPPS == ON || DUMPPS == ALL
	strcpy(filename, "psmat.ps");
	dump_ps_mat(filename, 0, 0, eval_size, eval_size, argv, argc, OPEN);
#endif
	//int mem_fma = memcount;
	mulMatDirect(sys); /* Compute the direct part matrices. */

#if DIRSOL == OFF		/* with DIRSOL just want to skip to solve */

#if PRECOND == BD
					   /* starttimer; */
	bdmulMatPrecond(sys);
	/* stoptimer */
	prsetime = dtime; /* preconditioner set up time */
#endif

#if PRECOND == OL
					  /* starttimer; */
	olmulMatPrecond(sys);
	/* stoptimer */
	prsetime = dtime; /* preconditioner set up time */
#endif

#if DMPREC == ON
	dump_preconditioner(sys, chglist, 1); /* dump prec. and P to matlab file */
#endif

#if DPSYSD == ON
	dissys(sys);
#endif

#if CKDLST == ON
	chkList(sys, DIRECT);
#endif

#endif				/* DIRSOL == OFF */
	dumpnums(ON, eval_size); /* save num/type of pot. coeff calcs */

	dirtimesav = dirtime; /* save direct matrix setup time */
	dirtime = 0.0; /* make way for direct solve time */

#if DIRSOL == OFF

#if DUMPPS == ON
	dump_ps_mat(filename, 0, 0, eval_size, eval_size, argv, argc, CLOSE);
#endif
	int mem_fma = memcount;

	/* starttimer; */
	//mulMatUp(sys); /* Compute the upward pass matrices. */
#if DNTYPE == NOSHFT
	mulMatDown(sys); /* find matrices for no L2L shift dwnwd pass */
#endif

#if DNTYPE == GRENGD
					 //mulMatDown(sys); /* find matrices for full Greengard dnwd pass*/
#endif

#if CKDLST == ON
	chkList(sys, DIRECT);
	chkLowLev(sys, DIRECT);
	chkEvalLstD(sys, DIRECT);
#endif

	//mulMatEval(sys); /* set up matrices for evaluation pass */ 


	mem_fma = memcount - mem_fma;
	/*FILE * ffpp;

	ffpp = fopen("file_mlaca_memory.txt", "w+");
	fprintf(ffpp, "  MLFMA memory: %d\n", mem_fma / 1024);
	fclose(ffpp);*/

	/* stoptimer */
	mulsetup = dtime; /* save multipole matrix setup time */

	dumpnums(OFF, eval_size); /* dump num/type of pot. coeff calcs */

#if DUMPPS == ALL
	dump_ps_mat(filename, 0, 0, eval_size, eval_size, argv, argc, CLOSE);
#endif

#if DISSYN == ON
	dumpSynop(sys);
#endif

#if DMTCNT == ON
	dumpMatBldCnts(sys);
#endif

#endif				/* DIRSOL == ON */
	/*
	fprintf(stdout, "\nITERATION DATA");
	ttliter = capsolve(&capmat, sys, chglist, eval_size, up_size, num_cond,
	name_list);
	*/
#if MKSDAT == ON		/* dump symmetrized, 4 pi eps scaled matrix */
	/*  mksCapDump(capmat, num_cond, relperm, &name_list); */
#endif

#if TIMDAT == ON 
	ttlsetup = initalltime + dirtimesav + mulsetup;
	multime = uptime + downtime + evaltime;
	ttlsolve = dirtime + multime + prectime + conjtime;

	if (indsys->opts->debug == ON)
		fprintf(stdout, "\nTIME AND MEMORY USAGE SYNOPSIS\n");
#endif

#ifdef OTHER
	if (TIMDAT == ON) {
		//fprintf(stdout,
		//"Warning: compilation with OTHER flag gives incorrect times\n");
	}
#endif

#if TIMDAT == ON
#if 1==0
	fprintf(stdout, "Total time: %g\n", ttlsetup + ttlsolve);
	fprintf(stdout, "  Total setup time: %g\n", ttlsetup);
	fprintf(stdout, "    Direct matrix setup time: %g\n", dirtimesav);
	fprintf(stdout, "    Multipole matrix setup time: %g\n", mulsetup);
	fprintf(stdout, "    Initial misc. allocation time: %g\n", initalltime);
	fprintf(stdout, "  Total iterative P*q = psi solve time: %g\n", ttlsolve);
	fprintf(stdout, "    P*q product time, direct part: %g\n", dirtime);
	fprintf(stdout, "    Total P*q time, multipole part: %g\n", multime);
	fprintf(stdout, "      Upward pass time: %g\n", uptime);
	fprintf(stdout, "      Downward pass time: %g\n", downtime);
	fprintf(stdout, "      Evaluation pass time: %g\n", evaltime);
	fprintf(stdout, "    Preconditioner solution time: %g\n", prectime);
	fprintf(stdout, "    Iterative loop overhead time: %g\n", conjtime);

	if (DIRSOL == ON) { /* if solution is done by Gaussian elim. */
		fprintf(stdout, "\nTotal direct, full matrix LU factor time: %g\n", lutime);
		fprintf(stdout, "Total direct, full matrix solve time: %g\n", fullsoltime);
		fprintf(stdout, "Total direct operations: %d\n", fulldirops);
	}
	else if (EXPGCR == ON) { /* if solution done iteratively w/o multis */
		fprintf(stdout, "\nTotal A*q operations: %d (%d/iter)\n",
			fullPqops, fullPqops / ttliter);
	}
#endif
	/*
	if (indsys->opts->debug == ON || true) {
	fprintf(stdout, "Total memory allocated: %ld kilobytes ",
	memcount / 1024);
	uallocEfcy(memcount);

	fprintf(stdout, "  Q2M  matrix memory allocated: %7.ld kilobytes\n",
	memQ2M / 1024);
	memcount = memQ2M;
	fprintf(stdout, "  Q2L  matrix memory allocated: %7.ld kilobytes\n",
	memQ2L / 1024);
	memcount += memQ2L;
	fprintf(stdout, "  Q2P  matrix memory allocated: %7.ld kilobytes\n",
	memQ2P / 1024);
	memcount += memQ2P;
	fprintf(stdout, "  L2L  matrix memory allocated: %7.ld kilobytes\n",
	memL2L / 1024);
	memcount += memL2L;
	fprintf(stdout, "  M2M  matrix memory allocated: %7.ld kilobytes\n",
	memM2M / 1024);
	memcount += memM2M;
	fprintf(stdout, "  M2L  matrix memory allocated: %7.ld kilobytes\n",
	memM2L / 1024);
	memcount += memM2L;
	fprintf(stdout, "  M2P  matrix memory allocated: %7.ld kilobytes\n",
	memM2P / 1024);
	memcount += memM2P;
	fprintf(stdout, "  L2P  matrix memory allocated: %7.ld kilobytes\n",
	memL2P / 1024);
	memcount += memL2P;
	fprintf(stdout, "  Q2PD matrix memory allocated: %7.ld kilobytes\n",
	memQ2PD / 1024);
	memcount += memQ2PD;
	fprintf(stdout, "  Miscellaneous mem. allocated: %7.ld kilobytes\n",
	memMSC / 1024);
	memcount += memMSC;
	fprintf(stdout, "  Inductance mem allocated: %7.ld kilobytes\n",
	memIND / 1024);
	memcount += memIND;
	fprintf(stdout, "  Total memory (check w/above): %7.ld kilobytes\n",
	memcount / 1024);
	}
	*/
#endif

	return sys;
}




