/***************************************************************************
  **************************************************************************
  
                           S2kit 1.0

          A lite version of Spherical Harmonic Transform Kit

   Peter Kostelec, Dan Rockmore
   {geelong,rockmore}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
   Copyright 2004 Peter Kostelec, Dan Rockmore

   This file is part of S2kit.

   S2kit is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   S2kit is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with S2kit; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/


/*****************************************************************

  FST_semi.c - routines common to FST_semi_fly.c and FST_semi_memo.c

  The primary functions in this package are

  4) TransMult() - Multiplies harmonic coefficients using Driscoll-Healy
                    result.  Dual of convolution in "time" domain.

  one utility function:

  1) seanindex(): Given bandwidth bw, seanindex(m,l,bw) will give the
     position of the coefficient f-hat(m,l) in the one-row array

  For descriptions on calling these functions, see the documentation
  preceding each function.


*/


#include <math.h>
#include <stdlib.h>
#include <s2kit/FST_semi.h>


/************************************************************************/


/*****************************************************************

  Given bandwidth bw, seanindex(m,l,bw) will give the position of the
  coefficient f-hat(m,l) in the one-row array that Sean stores the spherical
  coefficients. This is needed to help preserve the symmetry that the
  coefficients have: (l = degree, m = order, and abs(m) <= l)

  f-hat(l,-m) = (-1)^m * conjugate( f-hat(l,m) )

  Thanks for your help Mark!

  ******************************************************************/

int seanindex(int m,
	      int l,
	      int bw)
{
  int x = m * (m + 1) / 2;

  if( m >= 0 )
    return bw * m - x + l;
  return bw * (bw + m) + x + l;
}

/************************************************************************/
/*
  multiplies harmonic coefficients of a function and a filter.
  See convolution theorem of Driscoll and Healy for details.

  bw -> bandwidth of problem
  size = 2*bw

  datacoeffs should be output of an FST, filtercoeffs the 
  output of an FZT.  There should be (bw * bw) datacoeffs,
  and bw filtercoeffs.
  rres and ires should point to arrays of dimension bw * bw.

*/

void TransMult(double *rdatacoeffs, double *idatacoeffs, 
	       double *rfiltercoeffs, double *ifiltercoeffs, 
	       double *rres, double *ires,
	       int bw)
{

  int m, l, size;
  double *rdptr, *idptr, *rrptr, *irptr;

  size = 2*bw ;

  rdptr = rdatacoeffs;
  idptr = idatacoeffs;
  rrptr = rres;
  irptr = ires;

  for (m=0; m<bw; m++) {
    for (l=m; l<bw; l++) {
      compmult(rfiltercoeffs[l], ifiltercoeffs[l],
	       rdptr[l-m], idptr[l-m],
	       rrptr[l-m], irptr[l-m]);

      rrptr[l-m] *= sqrt(4*M_PI/(2*l+1));
      irptr[l-m] *= sqrt(4*M_PI/(2*l+1));

    }
    rdptr += bw-m; idptr += bw-m;
    rrptr += bw-m; irptr += bw-m;
  }
  for (m=bw+1; m<size; m++) {
    for (l=size-m; l<bw; l++){
      compmult(rfiltercoeffs[l], ifiltercoeffs[l],
	       rdptr[l-size+m], idptr[l-size+m],
	       rrptr[l-size+m], irptr[l-size+m]);

      rrptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));
      irptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));

    }
    rdptr += m-bw; idptr += m-bw;
    rrptr += m-bw; irptr += m-bw;
  }

}
