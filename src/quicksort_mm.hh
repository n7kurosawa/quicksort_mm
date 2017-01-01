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
// * thinning at the pivot selection[4].
//
// [1] C.A.R. Hoare, Commun. ACM 4, 321 (1961).
// [2] M. Blum, et al., J. Comput. Syst. Sci. 7, 448 (1973).
// [3] K. Chen, A. Dumitrescu, arXiv:1409.3600 [cs.DS] (2014).
// [4] N. Kurosawa, arXiv:1698.04852 [cs.DS] (2016).
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
  template<class RAIt>
  inline RAIt median3(RAIt a, RAIt b, RAIt c)
  {
    return (*a < *b)
      ? ((*b < *c) ? b : ((*a < *c) ? c : a))
      : ((*a < *c) ? a : ((*b < *c) ? c : b));
  }

  // Get the median of given five elements.
  template<class RAIt>
  inline RAIt median5(RAIt a, RAIt b, RAIt c, RAIt d, RAIt e)
  {
    if (*a > *b) std::swap(a, b);
    if (*c > *d) std::swap(c, d);
    if (*a > *c) {
      std::swap(a, c);
      std::swap(b, d);
    }
    if (*b > *e) std::swap(b, e);
    if (*b > *c) {
      return *b > *d ? d : b;
    }
    else {
      return *c > *e ? e : c;
    }
  }


  // Simple insertion sort
  template<class RAIt>
  inline void insertion_sort(RAIt first, RAIt last)
  {
    if (first == last) return;
    for (auto current = first+1; current != last; current++) {
      if (*current < *(current-1)) {
        auto cursor = current;
        auto target = *current;
        do {
          *cursor = *(cursor-1);
          cursor--;
        } while (first != cursor && target < *(cursor-1));
        *cursor = target;
      }
    }
  }


  // Hoare's Partition
  template<class RAIt>
  RAIt partition(RAIt first, RAIt last, RAIt pivot)
  {
    if (first != pivot) std::swap(*first, *pivot);
    pivot = first;
    auto lo = first, hi = last;
    for (;;) {
      for (;;) {
        hi--;
        if (lo == hi) goto PARTITION_END;
        if (*hi <= *pivot) break;
      }
      for (;;) {
        lo++;
        if (lo == hi) goto PARTITION_END;
        if (*lo >= *pivot) break; 
      }
      std::swap(*lo, *hi);
    }
  PARTITION_END:;
    std::swap(*pivot, *lo);
    return lo;
  }


  // ======================================================
  // Median of Medians
  // ======================================================

  template<int s_=2, class RAIt>
  RAIt rs3_5_2_select_pivot(RAIt first, RAIt last);

  template<int s_=2, class RAIt>
  RAIt rs3_5_2_select_kth(RAIt first, RAIt last, size_t k);


  // Extension of the repeated step algorithm (3-5).
  // 3-3 and 4-4 are presented in the original paper.
  template<int s_, class RAIt>
  RAIt rs3_5_2_select_pivot(RAIt first, RAIt last)
  {
    static_assert(s_ > 0, "Thinning parameter must be positive");
  
    size_t nelem = last - first;
    // The cutoff value 15 is taken as in the paper:
    // M. Durand, Inf. Process. Lett. 85, 73 (2003).
    if (nelem < 15) return first + nelem/2;
    // The following cutoff values are not optimized and should be refined.
    if (nelem < 80) return median3(first, first+nelem/2, last-1);
    if (nelem < std::max(30*s_,200)) return median5(first, first+nelem/4, first+nelem/2, first+3*nelem/4, last-1);

    size_t nnext = nelem/(15*s_);
    auto p = first + 0*(nelem/15);
    auto q = first + 7*(nelem/15);
    auto r = last - 7*nnext;

    for (size_t i = 0; i < nnext; i++) {
      auto x0 = median3(p+i*7+0, p+i*7+1, p+i*7+2);
      auto x1 = median3(p+i*7+3, p+i*7+4, p+i*7+5);
      auto x2 = median3(p+i*7+6, q+i,     r+i*7+0);
      auto x3 = median3(r+i*7+1, r+i*7+2, r+i*7+3);
      auto x4 = median3(r+i*7+4, r+i*7+5, r+i*7+6);

      auto xx = median5(x0, x1, x2, x3, x4);
      if (xx != q+i) std::swap(*xx, *(q+i));
    }
    return rs3_5_2_select_kth(q, q+nnext, nnext/2);
  }


  template<int s_, class RAIt>
  static RAIt rs3_5_2_select_kth(RAIt first, RAIt last, size_t k)
  {
    static_assert(s_ > 0, "Thinning parameter must be positive");

    size_t nelem = last - first;  
    if (nelem == 1) return first;
    if (nelem == 2) {
      if (first[0] > first[1]) std::swap(first[0], first[1]);
      return first+k;
    }

    auto pivot = rs3_5_2_select_pivot<s_>(first, last);
    auto pivotx = partition(first, last, pivot);
  
    size_t nl = pivotx - first;

    if (nl < k) {
      return rs3_5_2_select_kth(pivotx + 1, last, k-nl-1);
    }
    else if (k < nl) {
      return rs3_5_2_select_kth(first, first+nl, k);
    }
    else {
      return pivotx;
    }
  }


  // ======================================================
  // Quicksort with median of medians
  //
  // asymptotic comparison number for arrays of size N:
  // Random:  1.55 N ln N + O(N)
  // Worst:  21.33 N ln N + O(N)
  // ======================================================
  template<class RandomAccessIterator>
  void quicksort(RandomAccessIterator first, RandomAccessIterator last)
  {
    size_t nelem = last - first;

    if (nelem < 24) {
      insertion_sort(first, last);
      return;
    }

    auto pivot = rs3_5_2_select_pivot<21>(first, last);
    auto pivot_position = partition(first, last, pivot);

    if (last - pivot_position < pivot_position - first) {
      quicksort(pivot_position+1, last);
      quicksort(first, pivot_position);
    }
    else {
      quicksort(first, pivot_position);
      quicksort(pivot_position+1, last);
    }
  }


  // ======================================================
  // Quickselect with median of medians
  //
  // asymptotic comparison number for arrays of size N:
  // Random:  2.83 N + o(N)
  // Worst:  26.40 N + o(N)
  // ======================================================
  template<class RandomAccessIterator>
  void quickselect(RandomAccessIterator first, RandomAccessIterator kth, RandomAccessIterator last)
  {
    size_t k = kth - first;
    size_t nelem = last - first;
    if (nelem <= k) return;
    rs3_5_2_select_kth<21>(first, last, k);
  }
}


#endif
