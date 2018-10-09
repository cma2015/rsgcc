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

#ifndef MI_H
#define MI_H

#include "grid.h"
#include "points.h"

typedef struct {
  int k;
  int n;
  coord_t* psi;
  coord_t* sxs;
  int*     xiis;
  coord_t* sys;
  int*     yiis;
  grid_t grid;
} mi_t;


int make_mi(mi_t* const m, const int n, const int k);
void destroy_mi(mi_t* const m);
coord_t mutual_information(mi_t* const m, const coord_t* const xs, const coord_t* const ys);


coord_t c_gcc(mi_t* const m, const coord_t* const xs, const coord_t* const ys, const int* const xsix, const int* const ysix );
//coord_t gini_correlation(mi_t* const m, const coord_t* const xs, const coord_t* const ys);

coord_t c_pcc(mi_t* const m, const coord_t* const xs, const coord_t* const ys);
coord_t c_scc(mi_t* const m, const coord_t* const xs, const coord_t* const ys, const int* const xsix, const int* const ysix );
coord_t c_kcc(mi_t* const m, const coord_t* const xs, const coord_t* const ys );
coord_t c_eudist(mi_t* const m, const coord_t* const xs, const coord_t* const ys);


#endif //MI_H
