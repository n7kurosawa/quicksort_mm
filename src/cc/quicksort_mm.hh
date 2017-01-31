// Written in 2016 by n.kurosawa
//
// This program is under the CC0 Public Domain Dedication 1.0.
// See <http://creativecommons.org/publicdomain/zero/1.0/> for details. 
// This program is distributed without any warranty.


// ======================================================
// Quicksort/Quickselect with median of medians in C++ language.
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

#ifndef QUICKSORT_MM_HH_INCLUDED
#define QUICKSORT_MM_HH_INCLUDED

#include <algorithm>
#include <utility>

namespace quicksort_mm {
  // ======================================================
  // Utilities
  // ======================================================

  // Get the median of given three elements.
  template<class RAIt, class Cmp>
  inline RAIt median3(RAIt a, RAIt b, RAIt c, Cmp cmp)
  {
    return cmp(*a, *b)
      ? (cmp(*b, *c) ? b : (cmp(*a, *c) ? c : a))
      : (cmp(*a, *c) ? a : (cmp(*b, *c) ? c : b));
  }

  // Get the median of given five elements.
  template<class RAIt, class Cmp>
  inline RAIt median5(RAIt a, RAIt b, RAIt c, RAIt d, RAIt e, Cmp cmp)
  {
    if (cmp(*b, *a)) std::swap(a, b);
    if (cmp(*d, *c)) std::swap(c, d);
    if (cmp(*c, *a)) {
      std::swap(a, c);
      std::swap(b, d);
    }
    if (cmp(*e, *b)) std::swap(b, e);
    if (cmp(*c, *b)) {
      return cmp(*d, *b) ? d : b;
    }
    else {
      return cmp(*e, *c) ? e : c;
    }
  }


  // Simple insertion sort
  template<class RAIt, class Cmp>
  inline void insertion_sort(RAIt first, RAIt last, Cmp cmp)
  {
    if (first == last) return;
    for (auto current = first+1; current != last; current++) {
      if (cmp(*current, *(current-1))) {
        auto cursor = current;
        auto target = *current;
        do {
          *cursor = *(cursor-1);
          cursor--;
        } while (first != cursor && cmp(target, *(cursor-1)));
        *cursor = target;
      }
    }
  }


  // Hoare's Partition
  template<class RAIt, class Cmp>
  RAIt partition(RAIt first, RAIt last, RAIt pivot, Cmp cmp)
  {
    if (first != pivot) std::swap(*first, *pivot);
    pivot = first;
    auto lo = first, hi = last;
    for (;;) {
      for (;;) {
        hi--;
        if (lo == hi) goto PARTITION_END;
        if (!cmp(*pivot, *hi)) break;
      }
      for (;;) {
        lo++;
        if (lo == hi) goto PARTITION_END;
        if (!cmp(*lo, *pivot)) break; 
      }
      std::swap(*lo, *hi);
    }
  PARTITION_END:;
    std::swap(*pivot, *lo);
    return lo;
  }


  // approximate square root
  inline size_t approx_sqrt(size_t n)
  {
    int base = -1;
    while (0 < n) {
      base += 1;
      n /= 4;
    }
    return 1 << base;
  }


  // ======================================================
  // Median of Medians
  // ======================================================

  template<class RAIt, class Cmp>
  RAIt rs3_5_2_pick_pivot(RAIt first, RAIt last, Cmp cmp, size_t s=2);

  template<class RAIt, class Cmp>
  RAIt rs3_5_2_find_kth(RAIt first, RAIt last, size_t k, Cmp cmp, size_t s=2);


  // A variant of the repeated step algorithm (3-5).
  // 3-3 and 4-4 are presented in the original paper.
  template<class RAIt, class Cmp>
  RAIt rs3_5_2_pick_pivot(RAIt first, RAIt last, Cmp cmp, size_t s)
  {  
    if (s < 2) s = 2;
    size_t nelem = last - first;
    // The cutoff value 15 is taken as in the paper:
    // M. Durand, Inf. Process. Lett. 85, 73 (2003).
    if (nelem < 15) return first + nelem/2;
    // The following cutoff values are not optimized and should be refined.
    if (nelem < 80) return median3(first, first+nelem/2, last-1, cmp);
    if (nelem < s*30 || nelem < 200) return median5(first, first+nelem/4, first+nelem/2, first+3*nelem/4, last-1, cmp);

    size_t nnext = nelem/(15*s);
    auto p = first + 0*(nelem/15);
    auto q = first + 7*(nelem/15);
    auto r = last - 7*nnext;

    for (size_t i = 0; i < nnext; i++) {
      // median of 5 (median of 3)
      auto x0 = median3(p+i*7+0, p+i*7+1, p+i*7+2, cmp);
      auto x1 = median3(p+i*7+3, p+i*7+4, p+i*7+5, cmp);
      auto x2 = median3(p+i*7+6, q+i,     r+i*7+0, cmp);
      auto x3 = median3(r+i*7+1, r+i*7+2, r+i*7+3, cmp);
      auto x4 = median3(r+i*7+4, r+i*7+5, r+i*7+6, cmp);

      auto xx = median5(x0, x1, x2, x3, x4, cmp);
      if (xx != q+i) std::swap(*xx, *(q+i));
    }

    // Get the median of (pseudo-) medians 
    return rs3_5_2_find_kth(q, q+nnext, nnext/2, cmp, s);
  }


  template<class RAIt, class Cmp>
  static RAIt rs3_5_2_find_kth(RAIt first, RAIt last, size_t k, Cmp cmp, size_t s)
  {
    size_t nelem = last - first;

    if (nelem < 7) {
      insertion_sort(first, last, cmp);
      return first+k;
    }

    auto pivot = rs3_5_2_pick_pivot(first, last, cmp, s);
    auto pivotx = partition(first, last, pivot, cmp);
  
    size_t nl = pivotx - first;

    // Recursive application
    if (nl < k) {
      return rs3_5_2_find_kth(pivotx + 1, last, k-nl-1, cmp);
    }
    else if (k < nl) {
      return rs3_5_2_find_kth(first, first+nl, k, cmp);
    }
    else {
      return pivotx;
    }
  }


  // ======================================================
  // Quicksort with median of medians
  //
  // asymptotic comparison number for arrays of size N:
  // Random:  1.44 N ln N + o(N ln N)
  // Worst:  14.76 N ln N + o(N ln N)
  // ======================================================
  template<class RandomAccessIterator, class Compare>
  void quicksort_body(RandomAccessIterator first, RandomAccessIterator last, Compare cmp, size_t s)
  {
    size_t nelem = last - first;
    if (s < 10) s = 10;

    // Boundary condition
    if (nelem < 16) {
      insertion_sort(first, last, cmp);
      return;
    }

    // Partition
    auto pivot = rs3_5_2_pick_pivot(first, last, cmp, s);
    auto pivot_position = partition(first, last, pivot, cmp);

    // Recursive application
    // The tail call optimization is assumed
    // 12/17 ~ 0.7059 is an approximate value of sqrt(1/2)
    if (last - pivot_position < pivot_position - first) {
      quicksort_body(pivot_position+1, last, cmp, s*12/17);
      quicksort_body(first, pivot_position, cmp, s*12/17);
    }
    else {
      quicksort_body(first, pivot_position, cmp, s*12/17);
      quicksort_body(pivot_position+1, last, cmp, s*12/17);
    }
  }


  template<class RandomAccessIterator, class Compare>
  void quicksort(RandomAccessIterator first, RandomAccessIterator last, Compare cmp)
  {
    size_t nelem = last-first;
    quicksort_body(first, last, cmp, approx_sqrt(nelem));
  }


  template<class RandomAccessIterator>
  void quicksort(RandomAccessIterator first, RandomAccessIterator last)
  {
    std::less<typename std::iterator_traits<RandomAccessIterator>::value_type> cmp;
    quicksort(first, last, cmp);
  }


  // ======================================================
  // Quickselect with median of medians
  //
  // asymptotic comparision number for arrays of size N:
  // Random:  2.76 N + o(N)
  // Worst:  26.50 N + o(N)
  // ======================================================
  template<class RandomAccessIterator, class Compare>
  void quickselect(RandomAccessIterator first, RandomAccessIterator kth, RandomAccessIterator last, Compare cmp)
  {
    size_t k = kth - first;
    size_t nelem = last - first;
    if (nelem <= k) return;
    rs3_5_2_find_kth(first, last, k, cmp, approx_sqrt(nelem));
  }

  template<class RandomAccessIterator>
  void quickselect(RandomAccessIterator first, RandomAccessIterator kth, RandomAccessIterator last)
  {
    std::less<typename std::iterator_traits<RandomAccessIterator>::value_type> cmp;
    quickselect(first, kth, last, cmp);
  }
}


#endif
