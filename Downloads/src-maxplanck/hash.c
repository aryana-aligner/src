#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "hash.h"

//const int TABLESIZE=4087;

uint64_t uminus(uint64_t x, uint64_t y){
	if(x > y)
		return x - y;
	return y - x;
}

uint64_t hash_f1(uint64_t t)
{
	return  t%TABLESIZE;
	//return t;
}

uint64_t hash_f2(uint64_t t)
{
	return 1 + (t % (TABLESIZE - 1));
	//return 1;
}

uint64_t hash(uint64_t t,uint64_t attempt)
{
	//return (hash_f1(t)+attempt*hash_f2(t))%TABLESIZE;
	//return (t + attempt) & (TABLESIZE - 1);
	return (t + attempt) % (TABLESIZE);
}


void update_value(bwt_t *const bwt, uint64_t h_index,uint64_t value,uint64_t level, uint64_t index, int * best, int best_size, int * best_found, hash_element * table, uint64_t read_start, uint64_t read_size, uint64_t read_index, uint64_t groupid)
{
	if (table[h_index].level != level){
		table[h_index].value=0;
		//table[h_index].width = 0;
		table[h_index].index = 0;
	}
	//if(table[h_index].last == read_ind)
	//	return;
	table[h_index].groupid=groupid;
	//if(index < table[h_index].index)
	if(uminus(table[h_index].index, index) > 10 && table[h_index].value > 3 * value)
		return;
	table[h_index].index = index;
//	fprintf(stderr,"index :: %llu\n",index);
	table[h_index].value += value;
	//table[h_index].last = read_ind;


	//fprintf(stderr,"sssize %d\n",table[h_index].parts);
	table[h_index].match_start[table[h_index].parts]=read_start;
	table[h_index].match_index[table[h_index].parts]=read_index;
	table[h_index].matched[table[h_index].parts++]=read_size;
//	if (table[h_index].parts>5)
//		fprintf(stderr,"parts ::::::::::; %d\n",table[h_index].parts);
	//fprintf(stderr,"vuy %lld\n ",h_index);
//	fprintf(stderr," table :: %llu %llu\n",table[h_index].index,table[h_index].value);
//	fprintf(stderr,"before best updte %llu!\n", table[h_index].value);

	int i=0;
	int changed=0;
	for (i=0; i<best_size; i++)
		if (best[i]==(int)(h_index))
		{
			changed=1;
			int j=0;
			for (j=i+1; j<best_size; j++)
			{
				if (best[j]==-1 || table[best[j]].value < table[h_index].value)
					best[j-1]=best[j];
				else 
					break;
			}
			best[j-1]=(int)(h_index);
		}

		
	if (!changed)
		for (i=(*best_found); i<=best_size; i++)
		{
			if(i==best_size || (best[i]!=-1 && table[best[i]].value >= table[h_index].value) )
			{
				if (i!=0)
				{
					if ((*best_found)>0)// && i==(*best_found))
						(*best_found)--;
					best[i-1] = (int)(h_index);
				}
				break;
			}
			if (i>0)
				best[i-1]=best[i];
		}


	return;
	
	/*if(best->value <= table[h_index].value + (mate!=NULL ? mate->value : 0) )
	{
		
//		fprintf(stderr,"best updte %llu!\n",value);
		if (best->index!=table[h_index].index)		
		{
			if (pair_options->paired)
				(pair_options->second_mate)=(pair_options->best_mate);
			(*second)=(*best);
		}
		if (pair_options->paired)
			(pair_options->best_mate)=mate;
		(*best) = table[h_index];
	}
	else if (second->value<=table[h_index].value + (mate!=NULL ? mate->value : 0) )
	{
		(*second) = table[h_index];
	}
*/
	//fprintf(stderr, "inside update_value, best.index = %llu\n", best->index);
/*	if (read_start+read_size<table[h_index].match_start[table[h_index].parts-1])
	{
		table[h_index].match_start[table[h_index].parts]=read_start+read_size;
		table[h_index].match_index[table[h_index].parts]=read_index+read_size;
		table[h_index].matched[table[h_index].parts++]=0;
	}
	table[h_index].match_start[table[h_index].parts]=read_start;
	table[h_index].match_index[table[h_index].parts]=read_index;
	table[h_index].matched[table[h_index].parts++]=1;
*/
	//-------------------------------compute cigar
	//if (table[h_index].index > index)
//		table[h_index].cigar[
		
	//if(table[h_index].index - index < 10)
//	if(fabs(table[h_index].index - index) > 5)
//	if(read_ind != table[h_index].last)
}

void error(uint64_t level){
	fprintf(stderr,"Table Overflow. level = %"PRIu64"\n", level);
}
 
void add(bwt_t *const bwt, uint64_t place,uint64_t value,uint64_t level, uint64_t index, int * best, int best_size, int * best_found, hash_element * table, 
		uint64_t read_start, uint64_t read_size, uint64_t read_index, int whole_read_size,uint64_t groupid)
{
//	fprintf(stderr,"index :: %llu\n",index);
	uint64_t i=0;
	for (;i<TABLESIZE; i++)
	{
		uint64_t h_index=hash(place,i);
//		fprintf(stderr,"%d :: \n",h_index);
		if (table[h_index].level != level)
		{
			table[h_index].value = 0;
			table[h_index].place = place;
			//table[h_index].width = 0;
			table[h_index].level = level;
			table[h_index].index = index;
			table[h_index].last = -1;
			table[h_index].groupid = 0;

			table[h_index].parts = 0;
			table[h_index].mate = -1;
//			fprintf(stderr,"NEW ADD\n");
//			table[h_index].match_start[0]=whole_read_size;
		}
		if (table[h_index].level==level && table[h_index].place==place)
		{
		//	fprintf(stderr,"hhere!\n");
			if (table[h_index].groupid==groupid)
				break;
//			fprintf(stderr,"index :: %llu, window :: %llu h_index:: %llu, i :: %llu, groudpid :: %llu\n",index,place,h_index,i,groupid);
			update_value(bwt, h_index,value,level, index, best, best_size, best_found, table, read_start, read_size, read_index,groupid);
			break;
		}
	}
	if (i==TABLESIZE)
		error(level);
}



void reset_hash(hash_element  table[]){
	memset(table,-1,(sizeof (hash_element)) * (TABLESIZE));
	int i = 0;

	for (i=0; i<TABLESIZE; i++){
		table[i].place = table[i].value = table[i].level = table[i].index = 0;
		table[i].parts = 0;
		table[i].level=0;	
/*		table[i].match_start=(uint64_t *)malloc(300*(sizeof (uint64_t)));
		table[i].match_index=(uint64_t *)malloc(300*(sizeof (uint64_t)));
		table[i].matched=(uint64_t *)malloc(300*(sizeof (uint64_t)));*/
//		table[i].match_start=(uint64_t *)malloc(MAX_READ_SIZE*(sizeof (uint64_t)));
//		table[i].match_index=(uint64_t *)malloc(MAX_READ_SIZE*(sizeof (uint64_t)));
//		table[i].matched=(uint64_t *)malloc(MAX_READ_SIZE*(sizeof (uint64_t)));
	}
}

void print(hash_element * table)
{
	uint64_t i=0;
	for (i=0; i<TABLESIZE; i++){
		printf("%"PRIu64" , %d , %"PRIu64" \n",table[i].place,table[i].value,table[i].level);
		/*	if(table[i].value > maxvalue){
			maxvalue = table[i].value;
			index = table[i].place;
			}*/
	}
	//if(maxvalue > 0)
	//printf("index = %d, value = %d", index, maxvalue);
	printf("\n");
}


