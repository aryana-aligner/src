#include <stdio.h>
#include <inttypes.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <stdint.h>
#include <pthread.h>

#include "sam.h"
const int sam_line=5000;
char int_to_bp[4]={'A','C','G','T'};

int sam_generator(char *buffer, char *qname, int flag, char * rname, bwtint_t position, uint32_t mapq, char * cigar, char * rnext, bwtint_t pnext, long long tlen, ubyte_t *seq, ubyte_t * quality,int len)
{
	char  *seq_string=(char *) malloc(len+1),*quality_string=(char *) malloc(len+1);
	if(flag == 4)
		mapq = 0;
	seq_string[len]=quality_string[len]=0;
	int i;
	for (i=0; i<len; i++)
	{
		if (seq[i]<4)
			seq_string[i]=int_to_bp[seq[i]];
		else
			seq_string[i]='N';
		quality_string[i]=quality[i];
	}
	buffer[0] = '\0';
	uint64_t pos=position;
/*	if ((flag & 4) ==0)
		pos++;
	else 
		pos=0;
*/	int head=snprintf(buffer,sam_line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%"PRIu64"\t%lld\t%s\t%s\n",qname,flag,rname, pos,mapq, cigar,rnext,pnext,tlen,seq_string,quality_string);
	free (seq_string);
	free (quality_string);
	//fprintf(stderr, "%s\t%d\t%s\t%llu\t%lu\t%s\t%s\t%u\t%lld\t%s\t%s\n",qname,flag,rname,position + 1,mapq, cigar,"*",0,(long long)(0),seq_string,quality_string);
	return head;
}

int sam_headers(char * buffer, bwtint_t *  offset, int size){
	int head=0;
	int i;
	//fprintf(stderr, "offset[size] = %llu, offset[size - 1] = %llu, last = %llu, size = %d\n", offset[size], offset[size - 1], offset[size] - offset[size - 1], size);
	head+=snprintf(buffer,40,"@HD\tVN:1.0\tSO:unsorted\n");
	//fprintf(stderr, "buffer = %s offset[0] = %llu\n name[0] = %s\n", buffer, offset[0], name[0]);
	for(i=0; i<size; i++){ //TODO
		head+=snprintf(buffer+head,40,"@SQ\tSN:%s\tLN:%"PRIu64"\n",name[i], offset[i + 1] - offset[i]); //TODO
	}
	return head;
	//return buffer;
}

/*int main()
{
	char buf[400];
	//int head=0;
	char *names[2]={"ref1","ref2"};
	long long offset[2]={25,44};
	int head=sam_headers(buf,names,offset,2);
	sam_generator(buf,head,"qname",16,"ref",53,24,"cigar","rnext",999,19,"ACTF","#FDF#D");
	printf("%s",buf);
	return 0;
}
*/

