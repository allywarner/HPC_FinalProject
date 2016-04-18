#ifndef QUICKSORT_H_INCLUDED
#define QUICKSORT_H_INCLUDED
/*include guards */
#include <stdlib.h>

void quicksort (void* base, size_t num, size_t size,
                    int (*comp)(const void*, const void*));

#endif
