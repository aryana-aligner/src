//YA HAGh
#include "bwt.h"

#define max_refs 10000
//typedef int bwtint_t;
//const int sam_line=1000;

int load_reference(const char * filename);


int sam_generator(char *buffer, char *qname, int flag, char * rname, bwtint_t position, uint32_t mapq, char * cigar, char * rnext, bwtint_t pnext, long long tlen, ubyte_t *seq, ubyte_t * quality,int len);
char name[max_refs][1000];
//int sam_headers(char * buffer,char ** name, long long*  offset, int size);
int sam_headers(char * buffer, bwtint_t*  offset, int size);
