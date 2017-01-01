// Written in 2016-2017 by N.Kurosawa
//
// This program is under the CC0 Public Domain Dedication 1.0.
// See <http://creativecommons.org/publicdomain/zero/1.0/> for details. 
// This program is distributed without any warranty.

// Header of quicksort/quickselect with median of medians in C language.

#ifndef QUICKSORT_MM_H_INCLUDED_
#define QUICKSORT_MM_H_INCLUDED_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void quicksort_mm_quicksort(void *, size_t, size_t, int(const void *, const void *));
void quicksort_mm_quickselect(void *, size_t, size_t, size_t, int(const void *, const void *));

#ifdef __cplusplus
}
#endif

#endif
