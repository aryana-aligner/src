#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
//#include "bwtgap.h"
//#include "bwtaln.h"
#include "aligner.h"
//#include "utils.h"
//#include "bwt.h"
#include "smith.h"
#include <math.h>
#define false 0
#define true 1
#define bool int

const char atom[4]={'A','C','G','T'};

/*char getNuc(uint64_t place, uint64_t seq_len){
	int rev=0;
	if(place > (seq_len / 2))
	{
		place = (seq_len / 2) - (place - (seq_len / 2))-1;
		rev=1;
	}
	uint64_t block=place/(sizeof(bwtint_t)*4);
	int offset=place%(sizeof(bwtint_t)*4);
	uint64_t mask=3;
	mask=mask & (reference[block] >> (2*offset));
	if (rev==1)
		mask=3-mask;
	return mask;//atom[mask];
}*/

void showerr(int len, const ubyte_t *seq){
	int i = 0;
	for(i=0;i<len;i++)
		switch (seq[i]) {
			case 0:
				fprintf(stderr, "A");
				break;
			case 1:
				fprintf(stderr, "C");
				break;
			case 2:
				fprintf(stderr, "G");
				break;
			case 3:
				fprintf(stderr, "T");
				break;
			default:
				break;
		}
	fprintf(stderr, "\n");
}

void show(int len, const ubyte_t *seq){
	int i = 0;
	for(i=0;i<len;i++)
		switch (seq[i]) {
			case 0:
				printf("A");
				break;
			case 1:
				printf("C");
				break;
			case 2:
				printf("G");
				break;
			case 3:
				printf("T");
				break;
			default:
				break;
		}
	printf(" ");
}

void swap(ubyte_t *x, ubyte_t *y){
	ubyte_t tmp;
	tmp = *x;
	*x = *y;
	*y = tmp;
}
void swap64(uint64_t *x, uint64_t *y){
	uint64_t tmp;
	tmp = *x;
	*x = *y;
	*y = tmp;
}
void swapint(int *x, int *y){
	int tmp;
	tmp = *x;
	*x = *y;
	*y = tmp;
}


int create_cigar(hash_element *best, char *cigar, int len, const ubyte_t *seq, uint64_t seq_len,int **d, char **arr, char * tmp_cigar )
{
	int *valid=(int *)malloc(best->parts*(sizeof (int)));
	bwtint_t lastvalid=best->parts-1;
	bwtint_t i=0,j=0;
	for (i=best->parts-1; i>=0; i--)
	{
		valid[i]=1;
		if (abs((best->match_index[i]-best->match_start[i])-best->index) > 5+best->match_start[i]/20)
			valid[i]=0;
		if (valid[i]==1 && i<best->parts-1 && (best->match_start[lastvalid]+best->matched[lastvalid]>best->match_start[i]))
		{
			valid[i]=0;
			if (abs((best->match_index[i]-best->match_start[i])-best->index) < abs((best->match_index[lastvalid]-best->match_start[lastvalid])-best->index))
			{
				valid[lastvalid]=0;
				valid[i]=1;
			}
		}
		if (valid[i]==1)
			lastvalid=i;
		if (i==0)
			break;
	}

	//for (i=0; i<best->parts; i++)
	//	fprintf(stderr, "part %lld: %lld %lld %d %d\n",i,best->match_start[i],best->match_index[i],best->matched[i],valid[i]);

	int slack=20;
	bwtint_t head_match=0,head_index=best->index >= slack ? best->index-slack : 0;
	bwtint_t slack_index=head_index;
	int print_head=0;
	int total_errors=0;
	if (best->parts<=0 || best->parts > 50)
		fprintf(stderr , "too much parts!\n");
	if (0)
	{
		int errors_=0;
		print_head=smith_waterman(0, len, head_index, head_index+len+2*slack, cigar, print_head, seq, len, &errors_, seq_len, d, arr, tmp_cigar);
		total_errors+=errors_;
	}
	else
	for (i=best->parts-1; i>=0; i--)
	{
		if (valid[i] && !(best->match_start[i] < head_match || best->match_index[i] < head_index))// && abs((best->match_index[i]-head_index)-(best->match_start[i]-head_match)<=3+(best->match_index[i]-head_index)/20))
		{
			int errors=0;
			print_head=smith_waterman(head_match, best->match_start[i], head_index, best->match_index[i], cigar, print_head, seq, len, &errors, seq_len, d, arr, tmp_cigar);
			//fprintf(stderr,"start: %llu, end: %llu, errors :: %d\n",head_match,best->match_start[i],errors);
			total_errors+=errors;
			print_head+=snprintf(cigar+print_head,10,"%"PRIu64"%c",best->matched[i],'m');
			head_match=best->match_start[i]+best->matched[i];
			head_index=best->match_index[i]+best->matched[i];
		}
		if (i==0)
		{
			int errors=0;
			print_head=smith_waterman(head_match,len,head_index, head_index+len-head_match+slack, cigar, print_head,seq, len, &errors, seq_len, d, arr, tmp_cigar);
			//fprintf(stderr,"start: %llu, end: %llu, errors :: %d\n",head_match,len,errors);
			total_errors+=errors;
			break;
		}
	}

	//refining cigar
	int firstblood=1;
	int last_size=0;
	char last_char='m';
	print_head=0;
	for (i=0; cigar[i]; i++)
	{
		char tmp[10];
		j=0;
		while(isdigit(cigar[i]))
			tmp[j++]=cigar[i++];
		tmp[j]=0;
		bwtint_t num=atoi(tmp);
		if (firstblood==1)
		{
			firstblood=0;

			last_size=num;
			last_char=cigar[i];

			if (cigar[i]=='d')
			{
				last_size=0;
				slack_index+=num;
			}

			continue;
		}
		if (cigar[i]==last_char)
		{
			if (last_size>0)
				last_size+=num;
		}
		else
		{
			if (last_size>0)
				print_head+=snprintf(cigar+print_head,10,"%d%c",last_size,last_char);
			last_size=num;
			last_char=cigar[i];
		}
	}
	if (last_char!='d')
		print_head+=snprintf(cigar+print_head,10,"%d%c",last_size,last_char);

	best->index=slack_index;
	free(valid);
	//fprintf(stderr,"total errors:: %d\n",total_errors);
	return total_errors;//best->index;
}

void aligner(bwt_t *const bwt, int len, const ubyte_t *seq, bwtint_t level, hash_element * table, int *best, int best_size, int *best_found, aryana_args *args)
{

	//showerr(len,seq);
	//initialize
//	fprintf(stderr,"paired :: %d\n", options.paired);
	bwtint_t down, up;
	bwtint_t limit, i = 0, j = 0;
	
	/*best->index = best->value = best->place = best->level = 0;
	best->parts = 0;
	best->index = bwt->seq_len;
	second->index = second->value = second->place = second->level = 0;
	second->index = bwt->seq_len;
	second->parts=0;
*/
	//for (i=0; i<100; i++)
	//	fprintf(stderr,"%c",atom[getNuc(10000+i,bwt->seq_len)]);
	//fprintf(stderr,"\n");
	//set k
	bwtint_t k;
	if (args->seed_length==-1)
	{
		k = 26;
		if (len < 600)
			k=24;
		if (len < 300)
			k=22;
		if (len < 150)
			k=20;
		if (len <80)
			k=18;
		if (len <40)
		k=15;
	}
	else
		k=args->seed_length;

	//reversing
	j = len-1;
	i=0;
	for (; i<j && j >= 0 && i < len; (i++,j--)){
		swap(&seq[i],&seq[j]);
		//if (seq[i]>3 || seq[j]>3)
		//	return;
	}
	//inexact match
	bwtint_t groupid_last=1;
	for(i=len - 1;i>=k;i--){
		bwt_match_limit_rev(bwt, k, seq+i - k + 1, &down, &up,&limit);
		if(limit < k){
			i = i - k + limit;
			continue;
		}
		bwt_match_limit(bwt, i+1, seq, &down, &up,&limit);
		for (j=down; j<=up && j <= (down + 50); j++){
			bwtint_t index=bwt_sa(bwt,j);
			bwtint_t score=limit;

			bwtint_t rindex=index- (i - limit+1);
			add(bwt, rindex/len, score, level, index - (i - limit+1), best, best_size, best_found, table, i - limit+1, limit, index,len, groupid_last); // if level changed, check the find_value in hash.c
		}
		groupid_last++;
		if(i > k){
			if((limit - k + 1) > 0)
				i = i - limit + (k - 1);
			else
				fprintf(stderr, "manfi\n");
		}
	}
	if (args->best_factor!=-1)
	{
		for (i=best_size-2; i>=(*best_found); i--)
		{
			if (best[i] == -1)
				break;
			if (best[i]<best[best_size-1]/args->best_factor)
			{
				(*best_found)=i+1;
				break;
			}
			if (i==0)
				break;
		}
	}
		
}

