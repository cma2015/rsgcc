/*
Copyright 2010-2011 Gabriele Sales <gabriele.sales@unipd.it>


This file is part of parmigene.

knnmi is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License 
version 3 as published by the Free Software Foundation.

knnmi is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public
License along with parmigene. If not, see <http://www.gnu.org/licenses/>.
*/

#include "iqsort.h"
#include "mi.h"
#include <math.h>
#include <R.h>


static void init_psi(mi_t* const m) {
  static const coord_t c = 0.5772156649015328606065;

  m->psi = Calloc(m->n, coord_t);
  m->psi[0] = -c;

  int i;
  for (i = 1; i < m->n; i++)
    m->psi[i] = m->psi[i-1] + 1.0/i;
}

static coord_t get_psi(const mi_t* const m, const int i) {
  return m->psi[i-1];
}

#define compare_coords(a, b) ((**a) < (**b))

static void sort_coords(const coord_t* const cs, coord_t* const scs, int* const iis, const int n) {
  const coord_t* cps[n];

  int i;
  for (i = 0; i < n; i++)
    cps[i] = (coord_t*)(cs+i);

  QSORT(const coord_t*, cps, n, compare_coords);

  for (i = 0; i < n; i++) {
    scs[i] = *cps[i];
    iis[cps[i]-cs] = i;
  }
}

static dist_t find_range(const coord_t* const cs, const int center_idx, const int* const kis, const int k) {
  int i;
  dist_t md = 0;
  for (i = 0; i < k; i++) {
    const dist_t d = dist_abs(cs[center_idx] - cs[kis[i]]);
    if (d > md) md = d;
  }
  return md;
}

static int region_count(const coord_t* const scs, const int n, const int center_idx, const dist_t range) {
  const coord_t center = scs[center_idx];
  int c = 0;

  int i = center_idx-1;
  while (i >= 0 && center - scs[i] <= range) {
    c++;
    i--;
  }

  i = center_idx+1;
  while (i < n && scs[i] - center <= range) {
    c++;
    i++;
  }

  return c;
}

int make_mi(mi_t* const m, const int n, const int k) {
  if (n < k) return 0;

  m->k = k;
  m->n = n;
  init_psi(m);
  m->sxs  = Calloc(n, coord_t);
  m->xiis = Calloc(n, int);
  m->sys  = Calloc(n, coord_t);
  m->yiis = Calloc(n, int);

  return 1;
}

void destroy_mi(mi_t* const m) {
  Free(m->sxs);
  Free(m->xiis);
  Free(m->sys);
  Free(m->yiis);
  Free(m->psi);
}

coord_t mutual_information(mi_t* const m, const coord_t* const xs, const coord_t* const ys) {
  const coord_t* pxs;
  const coord_t* pys;
  make_grid(&m->grid, xs, ys, m->n, m->k);
  ordered_points(&m->grid, &pxs, &pys);

  sort_coords(pxs, m->sxs, m->xiis, m->n);
  sort_coords(pys, m->sys, m->yiis, m->n);

  int i;
  coord_t accum = 0;
  for (i = 0; i < m->n; i++) {
    int kis[m->k];
    search_knn(&m->grid, pxs[i], pys[i], kis);

    const dist_t mdx = find_range(pxs, i, kis, m->k);
    const int nx = region_count(m->sxs, m->n, m->xiis[i], mdx);

    const dist_t mdy = find_range(pys, i, kis, m->k);
    const int ny = region_count(m->sys, m->n, m->yiis[i], mdy);

    accum += get_psi(m, nx) + get_psi(m, ny);
  }

  destroy_grid(&m->grid);

  return get_psi(m, m->k) + get_psi(m, m->n) - (1.0/m->k) - (accum/m->n);
}


//for gini correlation
coord_t c_gcc(mi_t* const m, const coord_t* const xs, const coord_t* const ys, const int* const xsix, const int* const ysix ) {
  
  make_grid(&m->grid, xs, ys, m->n, m->k);
  
  const int len = m->n;
  int xorder, yorder;
  coord_t vecx_x[len], vecx_y[len], vecy_y[len], vecy_x[len];
  for (int i = 0; i < len; i++ ) {
  	xorder = xsix[i];
  	vecx_x[xorder-1] = xs[i];
  	vecy_x[xorder-1] = ys[i];
  	
  	yorder = ysix[i];
  	vecy_y[yorder-1] = ys[i];
  	vecx_y[yorder-1] = xs[i];
  }
  	

  coord_t accumx_x = 0;
  coord_t accumx_y = 0;
  coord_t accumy_y = 0;
  coord_t accumy_x = 0;
  for (int i = 0; i < len; i++) {
     //for order x
     coord_t weight = 2.0*(i+1)- m->n -1.0;
     accumx_x += weight*vecx_x[i];
     accumx_y += weight*vecx_y[i];
     
     accumy_y += weight*vecy_y[i];
     accumy_x += weight*vecy_x[i];
  }//for i

  coord_t gccx_y = accumx_y/accumx_x;
  coord_t gccy_x = accumy_x/accumy_y;
  coord_t final_gcc = gccx_y*gccx_y > gccy_x*gccy_x ? gccx_y : gccy_x;

  destroy_grid(&m->grid);

  return final_gcc;

}




//for euclidean distance
coord_t c_eudist( mi_t* const m, const coord_t* const xs, const coord_t* const ys) {
  int i;
  make_grid(&m->grid, xs, ys, m->n, m->k);

  coord_t sum1 = 0.0;
  coord_t t_tmp = 0.0;
  const int len = m->n;
  for( i = 0; i < len; i++ ) {
      t_tmp = xs[i] - ys[i];
      sum1 += t_tmp*t_tmp;
   }
    
  destroy_grid(&m->grid);

   if( sum1 == 0.0 ) return(0.0);
   else              return( sqrt(sum1) );
}




//for pearson correlation
coord_t c_pcc(mi_t* const m, const coord_t* const xs, const coord_t* const ys) {
  
  int i;
  make_grid(&m->grid, xs, ys, m->n, m->k);

  coord_t meanx = 0.0;
  coord_t meany = 0.0;
  coord_t sum1 = 0.0;
  coord_t sum2 = 0.0;
  coord_t sum3 = 0.0;
 
  const int len = m->n;
  for( i = 0; i < len; i++ ) {
      meanx += xs[i];
      meany += ys[i];
   }
   meanx = 1.0*meanx/len;
   meany = 1.0*meany/len;

   for( i = 0; i < len; i++ ) {
       sum1 += (xs[i] - meanx)*(ys[i] - meany);
       sum2 += (xs[i] - meanx)*(xs[i] - meanx);
       sum3 += (ys[i] - meany)*(ys[i] - meany);
    }
    
  destroy_grid(&m->grid);

   if( sum2 == 0.0 || sum3 == 0.0 ) return(0.0);
   else                             return( 1.0*sum1/(sqrt(sum2)*sqrt(sum3)) );

}

coord_t accsum( const int start, const int end) {
   return( 0.5*((end+1)*end - start*(start-1)) );
}


void maskrankforSCC( coord_t * valuevec, coord_t* ixvec, const int num ) {

  int i, j;
  int preIndex = 0;
  int lastIndex = 0;
  coord_t meanRank;
  for( i = 1; i < num; i++ ) {


     if( valuevec[i] != valuevec[i-1] ) {
        lastIndex = i - 1;
        
        if( preIndex < lastIndex ) {  //start to mask the rank
            meanRank = accsum( preIndex + 1, lastIndex + 1 )/(lastIndex - preIndex + 1);
            for( j = preIndex; j <= lastIndex; j++ ){
                  ixvec[j] = meanRank;
            }
        }
        //initialize 
        preIndex = i;
        lastIndex = 0;
 
     }//end if
  }//end for i


  if( preIndex >= 0 && preIndex < num - 1 ) {
     lastIndex = num - 1;
     meanRank = accsum( preIndex + 1, lastIndex + 1 )/(lastIndex - preIndex + 1);
     for( j = preIndex; j <= lastIndex; j++ )
        ixvec[j] = meanRank;
  }
}

//for spearman correlation
coord_t c_scc(mi_t* const m, const coord_t* const xs, const coord_t* const ys, const int* const xsix, const int* const ysix ) {
   make_grid(&m->grid, xs, ys, m->n, m->k);
  
  const int len = m->n;
  int i, j;
  
  //sort y by the rank of x
  int xorder, yorder;
  coord_t vecx_x[len], vecix_x[len], vecy_x[len], vecixy_x[len];
  for(i = 0; i < len; i++ ) {
     xorder = xsix[i];
     vecx_x[xorder-1] = xs[i];
     vecix_x[xorder-1] = xorder;
     vecy_x[xorder-1] = ys[i];
     vecixy_x[xorder-1] = ysix[i];
  }
  //mask the rank of x
  maskrankforSCC( vecx_x, vecix_x, len );
/*
  for( i = 0; i < len; i++ ) {
     printf("%f, %f, %f, %f\n", vecx_x[i], vecix_x[i], vecy_x[i], vecixy_x[i] );
  }
  printf("\n\n");
*/
  //sort x by the rank of y
  coord_t vecy_y[len], vecixy_y[len], vecixx_y[len];
  for( i = 0; i < len; i++ ) {
     yorder = vecixy_x[i];
     vecy_y[yorder-1] = vecy_x[i];
     vecixy_y[yorder-1] = yorder;
     vecixx_y[yorder-1] = vecix_x[i];
  }
  //mask the rank of y
  maskrankforSCC( vecy_y, vecixy_y, len );
/*
   for( i = 0; i < len; i++ ) {
     printf("%f, %f, %f\n", vecy_y[i], vecixy_y[i], vecixx_y[i] );
  }
  printf("\n\n");
*/
  coord_t tmp = 0;
  for( i = 0; i < len; i++ ) {
     tmp += (vecixy_y[i] - vecixx_y[i])*(vecixy_y[i] - vecixx_y[i]);
  }
//  printf("%d, %f\n", len, 1.0 - 6.0*tmp/(len*len*len - len) );
 
  destroy_grid(&m->grid);

  return( 1.0 - 6.0*tmp/(len*len*len - len) );
}



//for spearman correlation
coord_t c_kcc(mi_t* const m, const coord_t* const xs, const coord_t* const ys ) {
   make_grid(&m->grid, xs, ys, m->n, m->k);
  
  const int len = m->n;
  int i, j, num;
  num = 0;
  for( i = 1; i < len; i++ ) {
     for( j = 0; j < i; j++ ) {
         num += (xs[i] - xs[j])*(ys[i] - ys[j]) > 0 ? 1 : -1;
     }//end for j
  }//end for i

 
  destroy_grid(&m->grid);

  return( 2.0*num/(len*len - len) );
}


