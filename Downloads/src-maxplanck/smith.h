#include "bwt.h" 
//#include "bwa2.h"

int smith_waterman(uint64_t match_start, uint64_t match_end, uint64_t index_start, uint64_t index_end, char *cigar, int head, const ubyte_t *read, int len, int * errors, uint64_t seq_len
		,int **d, char **arr, char *tmp_cigar);
