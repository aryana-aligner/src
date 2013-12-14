#include "aryana_args.h"
#include "bwt.h"
//#include "bwtaln.h"
//#include "kbwt.h"
//#include "ttest.h"
#include "hash.h"

void aligner(bwt_t *const bwt, int len, const ubyte_t *seq, bwtint_t level, hash_element * table, int *best, int best_size, int *best_found, aryana_args *args);

int create_cigar(hash_element *best, char *cigar, int len, const ubyte_t *seq, uint64_t seq_len,int **d, char **arr, char * tmp_cigar );
//void aligner_test(bwt_t *const bwt, int len, const ubyte_t *seq, hash_element * table, uint64_t * kintervals);
//void aligner_old(bwt_t *const bwt, int len, const ubyte_t *seq, hash_element * table);

void showerr(int len, const ubyte_t *seq);

void show(int len, const ubyte_t *seq);
