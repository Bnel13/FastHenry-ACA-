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




typedef struct cx_struct {
	double real;
	double imag;
} CX;

#define cx_add(z,x,y) \
  do { \
    (z).real = (x).real + (y).real; \
    (z).imag = (x).imag + (y).imag; \
  } while(0)

#define cx_sub(z,x,y) \
  do { \
    (z).real = (x).real - (y).real; \
    (z).imag = (x).imag - (y).imag; \
  } while(0)

#define cx_mul(z,x,y) \
  do { \
    (z).real = (x).real * (y).real - (x).imag * (y).imag; \
    (z).imag = (x).imag * (y).real + (x).real * (y).imag; \
  } while(0)

#define cx_div(z,x,y) \
  do { \
    (z).real = 1.0 / ((y).real * (y).real + (y).imag * (y).imag); \
    (z).imag = (z).real; \
    (z).real *= (x).real * (y).real + (x).imag * (y).imag; \
    (z).imag *= (x).imag * (y).real - (x).real * (y).imag; \
  } while(0)

#define cx_conj_mul(z,x,y) \
  do { \
    (z).real = (x).real * (y).real + (x).imag * (y).imag; \
    (z).imag = (x).imag * (y).real - (x).real * (y).imag; \
  } while(0)

#define cx_abs(x) \
  (sqrt((x).real*(x).real+(x).imag*(x).imag))

#define cx_scalar_mult(z, alpha, x) \
  do { \
    (z).real = alpha*(x).real; \
    (z).imag = alpha*(x).imag; \
  } while(0)

static CX CXZERO = { 0, 0 };
static CX CXONE = { 1, 0 };
static CX CXMONE = { -1, 0 };

