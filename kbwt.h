#include <stdint.h>

#define MAX64 500000000

uint64_t before(uint64_t * kintervals, int size, uint64_t index);

uint64_t after(uint64_t * kintervals, int size, uint64_t index);

void next_interval(uint64_t * kintervals, int size, uint64_t *down, uint64_t *up);
