#include <stdio.h>

#include "kbwt.h"

inline int myflsll(uint64_t x) {
	/*int a = flsl(x >> 32);
	if (a) return a + 32;
	return flsl(x & 0xFFFFFFFF);*/
	return 0;
}

inline uint64_t before(uint64_t * kintervals, int size, uint64_t index){
	uint64_t i = index >> 6, s = size >> 6;
	int j = index & 63;
	if ((kintervals[i]>>j) & 1) return index;
	int f = myflsll(kintervals[i] & ((1 << (j + 1)) - 1));
	if (f) return ((i<<6) + f - 1);
	for (i--; i >=0; i--) {
		if(kintervals[i] == 0)
			continue;
		f = myflsll(kintervals[i]);
		if (f) return ((i<<6) + f - 1);
	}

	return 0;

	/*	//fprintf(stderr, "inside before => size = %d, index = %lld\n", size, index);


	//fprintf(stderr, "i = %lld, j = %lld, kintervals[%lld] = %x\n", i, j, i, kintervals[i]);
	int find = 0;
	while(!find){
		//fprintf(stderr, "i = %lld, j = %d, j<0: %d\n", i, j, (j < 0));
		if(j < 0){
			j = 63;
			i--;
		}
		if(i < 0)
			return 0;
		find = kintervals[i] >> (j) & 1;
		if(!find)
			j--;
	}
	//fprintf(stderr, "returning from before: i = %lld, j = %lld, result = %lld\n", i, j, ((i << 6) + j));
	return ((i << 6) + j);
	*/
}

inline uint64_t after(uint64_t * kintervals, int size, uint64_t index){
	uint64_t i = index >> 6, s = size >> 6;
	int j = index & 63;
	//if ((kintervals[i]>>j) & 1) return index;
	int f = ffsll(kintervals[i] >> (j + 1));
	if (f) return ((i<<6) + (j + 1) + f - 1);
	for (i++; i < s; i++) {
		if(kintervals[i] == 0)
			continue;
		f = ffsll(kintervals[i]);
		if (f) return ((i<<6) + f - 2);
	}

	return size - 1;
/*	int i = index >> 6;
	int j = index & 63;
	int find = 0;
	while(!find){
		if(j >= 64){
			j = 0;
			i++;
		}
		if(i > (size >> 6))
			return size - 1;
		if(((i << 6) + j) >= size)
			return -1;
		find = kintervals[i] >> (j) & 1;
		if(!find)
			j++;
	}
	if(find)
		j--;
	if(j < (index & 63) || i > (index >> 6))
		j = index & 63;
	return ((i << 6) + j);
	*/
}

void next_interval(uint64_t * kintervals, int size, uint64_t *down, uint64_t *up){
	//fprintf(stderr, "inside next_interval => size = %d, down = %lld, up = %lld\n", size, *down, *up);
	//if(*down > *up)
		//fprintf(stderr, "down > up\n");
	*down = before(kintervals, size, *down);
	//fprintf(stderr, "down completed. down = %lld\n", *down);
	*up = after(kintervals, size, *up);
	//fprintf(stderr, "up completed. up = %lld\n", *up);
}
