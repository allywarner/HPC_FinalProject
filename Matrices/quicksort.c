#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "quicksort.h"

void quicksort (void* base, size_t num, size_t size,
                  int (*comp)(const void*, const void*))
{
  if (num == 0 || num == 1)
    ;
  else {
  srand(time(NULL)); /* seed random number generator*/

  unsigned int pivot = rand() % num;
  unsigned int left, right;

  char temp[size];
  char* a = (char*)base;

  left = 0;
  right = num - 1;

  while (left < right) {
    if (left < pivot) {
      if (comp(a + left*size,a + pivot*size) <= 0)
        left++;
      else {
        memcpy(temp, a + left*size, size);
        memcpy(a + left*size, a + pivot*size, size);
        memcpy(a + pivot*size, temp, size);
        pivot = left;
      }
    } else {
      if (comp(a + right*size,a + pivot*size) >= 0)
        right--;
      else {
        memcpy(temp, a + right*size, size);
        memcpy(a + right*size, a + pivot*size, size);
        memcpy(a + pivot*size, temp, size);
        pivot = right;
      }
    }
  }

  quicksort(a, pivot, size, comp);
  quicksort(a + (pivot + 1)*size, num - (pivot + 1), size, comp);
}
}
