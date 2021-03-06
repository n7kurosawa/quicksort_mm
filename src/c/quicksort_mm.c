// Written in 2016-2017 by N.Kurosawa
//
// This program is under the CC0 Public Domain Dedication 1.0.
// See <http://creativecommons.org/publicdomain/zero/1.0/> for details. 
// This program is distributed without any warranty.

// ======================================================
// Quicksort/Quickselect with median of medians in C language.
//
// This routine combines the following ideas:
// * quicksort and quickselect[1], of course.
// * median of medians or BFPRT[2].
// * repeated step algorithms[3].
// * square root sampling[4].
// * thinning at the pivot selection[5].
//
// [1] C.A.R. Hoare, Commun. ACM 4, 321 (1961).
// [2] M. Blum, et al., J. Comput. Syst. Sci. 7, 448 (1973).
// [3] K. Chen, A. Dumitrescu, arXiv:1409.3600 [cs.DS] (2014).
// [4] C.C. McGeoch, J.D. Tygar, Random Struct. Alg. 7, 287 (1995).
// [5] N. Kurosawa, arXiv:1698.04852 [cs.DS] (2016).
// ======================================================



#include <assert.h>
#include <stdint.h>
#include "quicksort_mm.h"

typedef int (*comparator)(const void *, const void *);


// ======================================================
// Utilities
// ======================================================


// Swap the contents of the pointers.
// This routine is not efficient. We take simplicity.
static void swap(char *p, char *q, size_t sz)
{
  while (sz > 0) {
    char tmp = *p;
    *p = *q;
    *q = tmp;
    p++;
    q++;
    sz--;
  }
}


// Swap the contents of the pointers, unless the pointers are same. 
static void swap_unless_same(char *p, char *q, size_t sz)
{
  if (p != q) swap(p, q, sz);
}


// Swap two void pointer variables.
static void swap_ptr(char **p, char **q)
{
  char *tmp = *p;
  *p = *q;
  *q = tmp;
}




// Get the median of given five elements.
static char *median5(char *a, char *b, char *c, char *d, char *e, comparator cmp)
{
  if (cmp(a, b) > 0) swap_ptr(&a, &b);
  if (cmp(c, d) > 0) swap_ptr(&c, &d);

  if (cmp(a, c) > 0) {
    swap_ptr(&a, &c);
    swap_ptr(&b, &d);
  }
  if (cmp(b, e) > 0) swap_ptr(&b, &e);

  if (cmp(b, c) > 0) {
    return (cmp(b, d) > 0 ? d : b);
  }
  else {
    return (cmp(c, e) > 0 ? e : c);
  }
}


// Get the median of given three elements.
static inline char *median3(char *p, char *q, char *r, comparator cmp)
{
  return cmp(p, q) < 0
    ? (cmp(q, r) < 0 ? q : (cmp(p, r) < 0 ? r : p))
    : (cmp(p, r) < 0 ? p : (cmp(q, r) < 0 ? r : q));
}


// approximate square root
static inline size_t approx_sqrt(size_t n)
{
  int base = -1;
  while (0 < n) {
    base += 1;
    n /= 4;
  }
  return 1 << base;
}





// Hoare's Partition
static char *partition(char *begin, char *pivot, size_t n, size_t sz, comparator cmp)
{
  // Cursors
  char *lo = begin;
  char *hi = begin + sz * n;

  assert((hi-lo) % sz == 0);
  assert((pivot-begin) % sz == 0);
  assert((pivot-begin) / sz < n);

  // Setup pivot
  // For simplicity, we place pivot element at the first.
  // (limited version of 3-way partition).
  swap_unless_same(pivot, lo, sz);
  pivot = lo;

  // Partition
  for (;;) {
    int cmp_hi, cmp_lo;
    for (;;) {
      hi-=sz;
      if (lo >= hi) goto PARTITION_END;
      cmp_hi = cmp(hi, pivot);
      if (cmp_hi <= 0) break;
    }
    for (;;) {
      lo+=sz;
      if (lo >= hi) goto PARTITION_END;
      cmp_lo = cmp(lo, pivot);
      if (cmp_lo >= 0) break;
    }

    if (cmp_lo != 0 || cmp_hi != 0) swap(lo, hi, sz);
  }
PARTITION_END:;
  swap(pivot, lo, sz);
    
  assert(lo-begin >= 0);
  //assert(end-lo >= 0);

  return lo;
}



// ======================================================
// Median of Medians
// ======================================================

static char *rs3_5_2_pick_pivot(char *p, size_t n, size_t sz, size_t thin, comparator cmp);
static char *rs3_5_2_find_kth(char *p, size_t n, size_t sz, size_t thin, size_t kth, comparator cmp);

// A variant of the repeated step algorithm (3-5).
// 3-3 and 4-4 are presented in the original paper.
static char *rs3_5_2_pick_pivot(char *p, size_t n, size_t sz, size_t thin, comparator cmp)
{
  if (thin < 2) thin = 2;
  if (n < 15) return p + (n/2)*sz;
  if (n < 80) return median3(p, p+(n/2)*sz, p+(n-1)*sz, cmp);
  if (n < 30*thin || n < 200) return median5(p, p+(n/4)*sz, p+(n/2)*sz, p+(3*n/4)*sz, p+(n-1)*sz, cmp);    

  size_t nnext = n/(15*thin);
  char *p0 = p;
  char *q0 = p + 7*(n/15)*sz;
  char *r0 = p + (n-nnext*7)*sz;

  for (size_t i = 0; i < nnext; i++) {
    // median of 5 (median of 3)
    char *s0 = median3(p0+(i*7+0)*sz, p0+(i*7+1)*sz, p0+(i*7+2)*sz, cmp);
    char *s1 = median3(p0+(i*7+3)*sz, p0+(i*7+4)*sz, p0+(i*7+5)*sz, cmp);
    char *s2 = median3(p0+(i*7+6)*sz, q0+(i*1+0)*sz, r0+(i*7+0)*sz, cmp);
    char *s3 = median3(r0+(i*7+1)*sz, r0+(i*7+2)*sz, r0+(i*7+3)*sz, cmp);
    char *s4 = median3(r0+(i*7+4)*sz, r0+(i*7+5)*sz, r0+(i*7+6)*sz, cmp);

    swap_unless_same(q0+i*sz, median5(s0,s1,s2,s3,s4,cmp), sz);
  }

  // Get the median of (pseudo-) medians 
  return rs3_5_2_find_kth(q0, nnext, sz, 2, nnext/2, cmp);
}


static char *rs3_5_2_find_kth(char *p, size_t n, size_t sz, size_t thin, size_t kth, comparator cmp)
{
  assert(kth < n);

  if (n == 1) return p;
  if (n == 2) {
    if (cmp(p, p+sz) > 0) swap(p, p+sz, sz);
    return p+kth*sz;
  }

  char *pivot = rs3_5_2_pick_pivot(p, n, sz, thin, cmp);

  assert(p <= pivot);
  assert((pivot - p) % sz == 0);
  assert((pivot - p) / sz < n);
  
  char *pivotx = partition(p, pivot, n, sz, cmp);

  assert((pivotx - p) % sz == 0);
  assert((pivotx - p) / sz < n);
  assert(pivotx >= p);

  size_t nl = (pivotx - p) / sz;
  assert(n >= nl+1);
  size_t nr = n - nl - 1;

  // Recursive application
  if (nl < kth) {
    return rs3_5_2_find_kth(pivotx + sz, nr, sz, 2, kth-nl-1, cmp);
  }
  else if (kth < nl) {
    return rs3_5_2_find_kth(p, nl, sz, 2, kth, cmp);
  }
  else {
    return pivotx;
  }
}






// ======================================================
// Main routine of quick sort
// ======================================================

static void quicksort_body(char *begin, char *end, size_t sz, comparator cmp, size_t thin)
{
  assert(begin <= end);
  assert(sz > 0);
  assert((size_t)(end - begin) % sz == 0);

  if (thin < 10) thin = 10;

  // Length of the input array
  size_t n = (size_t)(end - begin) / sz;

  // Boundary condition
  if (n <= 1) {
    return;
  }
  if (n == 2) {
    if (cmp(begin, begin+sz) > 0) {
      swap(begin, begin+sz, sz);
    }
    return;
  }

  // Partition
  char *pivot = rs3_5_2_pick_pivot(begin, n, sz, thin, cmp);
  char *pivot_pos = partition(begin, pivot, n, sz, cmp);

  // Recursive application
  // The tail call optimization is assumed
  // 12/17 ~ 0.7059 is an approximate value of sqrt(1/2)
  if (end-pivot_pos < pivot_pos-begin) {
    quicksort_body(pivot_pos+sz, end, sz, cmp, thin*12/17);
    quicksort_body(begin, pivot_pos, sz, cmp, thin*12/17);
  }
  else {
    quicksort_body(begin, pivot_pos, sz, cmp, thin*12/17);
    quicksort_body(pivot_pos+sz, end, sz, cmp, thin*12/17);
  }
}



// ======================================================
// Quicksort with median of medians
//
// asymptotic comparison number for arrays of size N:
// Random:  1.44 N ln N + o(N ln N)
// Worst:  14.76 N ln N + o(N ln N)
// ======================================================
void quicksort_mm_quicksort(void *p, size_t n, size_t sz, comparator cmp)
{
  if (!p) return;
  if (n == 0) return;
  if (sz == 0) return;
    
  char *begin = (char *)p;
  char *end = begin + n*sz;

  if (end < begin) return; // In this case the routine does not work.
  quicksort_body(begin, end, sz, cmp, approx_sqrt(n));
}



// ======================================================
// Quickselect with median of medians
//
// asymptotic comparison number for arrays of size N:
// Random:  2.76 N + o(N)
// Worst:  26.50 N + o(N)
// ======================================================
void quicksort_mm_quickselect(void *p, size_t n, size_t sz, size_t kth, comparator cmp)
{
  if (!p) return;
  if (sz == 0) return;
  if (n == 0) return;
  if (n <= kth) return;
    
  char *begin = (char *)p;
  char *end = begin + n*sz;

  if (end < begin) return; // In this case the routine does not work.
  rs3_5_2_find_kth(p, n, sz, approx_sqrt(n), kth, cmp);
}
