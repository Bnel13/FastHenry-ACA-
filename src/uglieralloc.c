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



/* # ***** sort to /src/misc
 # ***** */

/* A memory allocator with only the ability to free all memory allocated.
 It is a modification of the description below.    MK 11/95
 */

/* 
 memory allocator for fastcap
 - almost identical to Kerigan & Ritchie sec 8.7
 - differs in that there is no free list since fastcap never frees any memory
 - also the amount of memory sbrk()'d is exactly equal (within a unit) to 
 the requested amount if its more than NALLOC units
 this cuts down on unused blocks usually added to the free list
 while still keeping the number of sbrk()'s down
 - the regular allocation functions are still availible for stdio to use
 although  stdio buffers are passed over when new memory is set up
 - machine dependancy should be confined to the type `ALIGN' and
 `MORECORE()' (i.e. non-UNIX machines won't have brk(), sbrk())
 - no attempt is made to make allocation efficient in terms of virtual pages
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef OTHER
#include <unistd.h>
#endif // OTHER

/* SRW */
#define NO_SBRK

#define NALLOC 8184		/* >= sizeof(HEADER)*NALLOC bytes sbrk()'d */
#define MAGICN 0xaaaaaaaaL	/* used to check fidelity of allocated blks */
#define UGDEBG 0		/* 0=> most efficient, no verification
				   1=> ualloc_verify() checks magicno, size
				   2=> ualloc_verify() prints lots to stdout
				   - can reduce efficiency due to big HEADER */

typedef double ALIGN; /* all allocated mem starts on an ALIGN bdry */

union header { /* allocated block header - all mem allocated
 in units of sizeof(HEADER) bytes */
	struct {
#if UGDEBG == 1 || UGDEBG == 2
		long magicno; /* to check for overwrites */
		union header *ptr; /* next allocated block */
#endif
#if UGDEBG == 2
		union header *enlarged; /* limit of newly enlarged block */
		union header *enlgdfrom; /* base pointer of newly added block */
		unsigned int osize; /* size before enlargement/truncation */
		unsigned int request; /* nbytes requested from allocator */
#endif
		unsigned int size; /* block size (memory + 1 or 2 header units) */
	} s;
	union header *lastblock; /* For freeing (MK  11/95) */
	ALIGN x; /* for alignment only, never referenced */
};

typedef union header HEADER;

/* 
 uses calloc to get a new block of memory
 - any header space put in by calloc is wasted
 - zeros memory
 - an alternative to mocore() but should only be used if sbrk() doesnt zero
 */
#define MORECORE(SIZE) (HEADER *)calloc(1, SIZE*sizeof(HEADER))

static HEADER *base = NULL; /* base of allocated block list */
static HEADER *allocp = NULL; /* last allocated block */
static unsigned int sizeofHDR = sizeof(HEADER);

static HEADER *lastblock = NULL; /* pointer to last allocated block */
/* for freeing purposes */

/* SRW */
static HEADER *mocore(unsigned int);
void fill_it(void*, char, int);
void ufree(void);
char *ualloc(unsigned int);
void ualloc_verify(void);
void uallocEfcy(long);

/*
 asks operating system for more memory which is added to the top block
 - memory not zeroed out
 */
static HEADER *mocore(unsigned int nu) {
	HEADER *cp;

	cp = (HEADER *) calloc(nu + 1, sizeofHDR);
	if (cp == NULL)
		return (NULL);

	/* set pointer to last block for freeing */
	cp->lastblock = lastblock;
	lastblock = cp;

#ifdef MATTDEBUG2
	printf("%x\n", cp);
#endif

	return (cp + 1);
}

/* fills space with value */
void fill_it(void *mem, char k, int nbytes) {
	/* printf("will fill: %x + 1 = %x\n",(HEADER *)mem - 1, (HEADER *)mem);*/
	memset(mem, k, nbytes);
}

/* frees lastblock and  the linked list of each block before it
 */
void ufree(void) {
	HEADER *ptr, *next;

	ptr = lastblock;
	while (ptr != NULL) {
		next = ptr->lastblock;

#ifdef MATTDEBUG2
		printf("%x\n", ptr);
#endif

		free(ptr);
		ptr = next;
	}

	allocp = lastblock = NULL;

}

/* 
 ugly storage allocator
 - no frees
 - allocates in blocks at least NALLOC*sizeof(HEADER) bytes long
 - ultimately uses mocore(), since no frees are done (sbrk() zeros added
 memory) this allocator performs like calloc() w/no explicit assigns to 0
 */
char *ualloc(unsigned int nbytes) {
	HEADER *p, *q;
	int nunits; /* size in number of sizeof(HEADER)'s */
	int brkunits; /* number to add to heap */
#if UGDEBG == 2
	HEADER *dummy = 0x0;
	HEADER *dummyret = 0x0;
#endif
	static unsigned int topblksize; /* size of current top block */

	if (nbytes == 0)
		return (NULL); /* should probably return something else */

#if UGDEBG == 1 || UGDEBG == 2
	nunits = 3+(nbytes-1)/sizeofHDR; /* rm for 2 hdrs too */
#else
	nunits = 2 + (nbytes - 1) / sizeofHDR; /* rm for 1 hdr (gets subtracted later)*/
#endif

	if ((q = allocp) == NULL) { /* no allocation yet */
		if ((q = allocp = base = mocore(1)) == NULL)
			return (NULL);
#ifdef MATTDEBUG
		printf("started allocating:  size of header == %d\n",sizeof(HEADER));
#endif
		topblksize = 1;
#if UGDEBG == 1 || UGDEBG == 2
		base->s.size = 1;
		base->s.ptr = base;
		base->s.magicno = MAGICN;
#endif
	}

	/*
	 check previously allocated block for room
	 - if it's big enough, use head (not tail) part of it
	 - if it's not, break out more memory and discard the previous block
	 */
	if (topblksize >= nunits) { /* if it's big enough */
#if UGDEBG == 1 || UGDEBG == 2
		p = q; /* copy old top block pointer */
		q += (nunits - 1); /* make q new top block pointer */
		q->s.size = p->s.size - nunits + 1; /* upper block size */
#else
		allocp = q + nunits - 1;
#endif
		topblksize -= (nunits - 1);
#if UGDEBG == 2
		p->s.osize = p->s.size; /* save size before new block carved out */
		p->s.request = nbytes; /* actual #bytes requested */
		p->s.enlarged = dummy;
		p->s.enlgdfrom = dummyret;
#endif
#if UGDEBG == 1 || UGDEBG == 2
		q->s.ptr = base; /* top block connected back to base */
		p->s.ptr = q; /* old top block pointer to new top block */
		if(allocp != p) allocp->s.ptr = p; /* link previous old block */
		q->s.magicno = MAGICN;
		p->s.size = nunits; /* set up returned block size */
		allocp = q; /* save the top block pointer */
		return((char *) (p+1)); /* return pntr to first unit beyond header */
#else
		return ((char *) q); /* use old header location */
#endif
	} else { /* get more memory, add to top block beyond
	 any stdio buffers made since last ualloc */
		brkunits = (nunits > NALLOC) ? nunits : NALLOC;
		if ((q = mocore(brkunits)) == NULL)
			return (NULL);
#if UGDEBG == 1 || UGDEBG == 2
		q->s.size = brkunits; /* setup size, etc. of new top block */
		q->s.magicno = MAGICN;
#endif
		topblksize = brkunits;
#if UGDEBG == 2
		dummy = q + q->s.size - 1; /* limit pointer of new block */
		dummyret = q; /* base pointer of new block */
#endif
		/* repeat above code to aviod having a loop */
#if UGDEBG == 1 || UGDEBG == 2
		p = q;
		q += (nunits - 1); /* make q new top block pointer */
		q->s.size = p->s.size - nunits + 1; /* upper block size */
#else
		allocp = q + nunits - 1;
#endif
		topblksize -= (nunits - 1);
#if UGDEBG == 2
		p->s.osize = p->s.size; /* save size before new block carved out */
		p->s.request = nbytes; /* actual #bytes requested */
		p->s.enlarged = dummy;
		p->s.enlgdfrom = dummyret;
#endif
#if UGDEBG == 1 || UGDEBG == 2
		q->s.ptr = base; /* top block connected back to base */
		p->s.ptr = q; /* old top block pointer to new top block */
		if(allocp != p) allocp->s.ptr = p; /* link previous old block */
		q->s.magicno = MAGICN;
		p->s.size = nunits; /* set up returned block size */
		allocp = q; /* save the top block pointer */
		return((char *) (p+1)); /* return pntr to first unit beyond header */
#else
		return ((char *) q); /* use old header location */
#endif
	}
}

/*
 checks for overwrites at end of blocks
 - checks magic numbers and sizes (UGDEBG == 1 or 2)
 - checks if length corresponds to pointers (UGDEBG == 1 or 2)
 - prints information about list of allocated blocks (UGDEBG == 2)
 */
void ualloc_verify(void) {
	HEADER *p;
	int cnt = 1;

#if UGDEBG == 1 || UGDEBG == 2
	for(p = base;; p = p->s.ptr, cnt++) {
		if(p == base && cnt > 1) break;
#endif
#if UGDEBG == 2
	fprintf(stdout,
			"%d 0x%x 0x%x %u %u bytes (osize %u enlarged from 0x%x to 0x%x)\n",
			cnt, p, p->s.ptr, p->s.size, p->s.request,
			p->s.osize, p->s.enlgdfrom, p->s.enlarged);
#endif
#if UGDEBG == 1 || UGDEBG == 2
	if(p->s.magicno != MAGICN) {
		fflush(stdout);
		fprintf(stderr, "ualloc_verify: bad block %d magic number\n", cnt);
		fflush(stderr);
		abort();
	}
	if(p->s.ptr - p + 1 < p->s.size) {
		fflush(stdout);
		fprintf(stderr, "ualloc_verify: bad block size ", cnt);
		fprintf(stderr, "(from 0x%x to 0x%x size %u)\n", p, p->s.ptr, p->s.size);
		fflush(stderr);
		abort();
	}
}
#endif
}

/*
 compares the total address range broken out to amount requested
 - efficiency is 100*memcount/(sbrk(0)-base) = 100*requested/(ttl broken out)
 - UGDEBG == 2 prints waste values:
 ualloc waste = 100*(mem discarded)/(ttl broken out)
 stdio waste = 100*(mem discarded to apease stdio usage)/(ttl broken out)
 (stdio waste is often zero and ualloc waste does not include the
 memory lost in each header struct (a problem with many small things)
 - if base == NULL (not using ugly allocator), final break value is printed
 */
void uallocEfcy(long memcount) {
#ifdef NO_SBRK
	return;
#else
#if UGDEBG == 2
	HEADER *p;
	unsigned int waste = 0;
	unsigned int stdiowaste = 0;
	int first = 1;
#endif
	int total;

	total = (int)((char*)sbrk(0) - (char*)base);

	if(base == NULL) fprintf(stdout, "(top of memory = 0x%lx", (long)sbrk(0));
	else fprintf(stdout, "(%.3g%% efficiency",
			100*((double)memcount)/((double)total));

#if UGDEBG == 2
	/* loop through blocks, get total btyes wasted (by allocator) */
	for(p = base;; p = p->s.ptr) {
		if(p == base && first == 0) break;
		if(p == NULL) break; /* compatability when ualloc not used */
		if(p->s.request == 0) {
			if((p->s.ptr)->s.size >= p->s.size) { /* if block dropped bec. too sm. */
				if(p != allocp) waste += p->s.size - 1;
				else waste += p->s.size;
			}
			else {
				if(p != allocp) stdiowaste += p->s.size - 1;
			}
		}
		first = 0;
	}
	if(p != NULL)
	fprintf(stdout,
			", waste: %.2g%% (ualloc), %.2g%% (stdio)",
			100*(double)(waste*sizeof(HEADER))/((double)total),
			100*(double)(stdiowaste*sizeof(HEADER))/((double)total));
#endif
	fprintf(stdout, ")\n");
#endif
}
