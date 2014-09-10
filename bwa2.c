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
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>

#include "bntseq.h"
#include "utils.h"
//#include "main.h"
#include "bwt.h"

#include "bwa2.h"
//#include "aligner.h"
#include "sam.h"

//#define num_of_threads 6
//#define read_size 0x44000
#define read_size 600000
#define each_read_size 5000
#define max_sam_line 3*MAX_READ_SIZE+1000
#define MAX_CIGAR_SIZE 6000
#define output_buffer (1 << 9)
#define sleep_time (1)
#define true 1
#define false 0

int counter = 0;
double base=10.0;
double scmax=100.0;
bwa_seq_t *seqs;
int n_seqs = 0;
int tot_seqs = 0;
bwtint_t offset[1000];
bwtint_t offInd = 0;
bwt_t *bwt;
//int head[num_of_threads];
//int tail[num_of_threads];
pthread_mutex_t input = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t output = PTHREAD_MUTEX_INITIALIZER;
//char output[num_of_threads][output_buffer][max_sam_line];
bwa_seqio_t *ks;
bwa_seqio_t *ks1,*ks2;
//char out_buf[num_of_threads+20][2*max_sam_line];
//int finish[num_of_threads];
//#define is_finished (finish[0] && finish[1] && finish[2] && finish[3] && finish[4] && finish[5] && finish[6] && finish[7] && finish[8] && finish[9] && finish[10] && finish[11])
uint32_t reference_chars;
uint64_t * reference;
uint32_t  reference_size;
uint32_t  reference_reminder;
uint64_t * ref;
uint32_t refSize;

int offsearch(int bot, int top, uint64_t * a, uint64_t key){
	if(bot > top)
		return -1;
	int mid = (bot + top) / 2;
	//fprintf(stderr, "bot = %d, top = %d, a[%d] = %llu\n", bot, top, mid, a[mid]);
	if(key >= a[mid] && a[mid + 1] > key)
		return mid;
	if(key < a[mid])
		return offsearch(bot, mid - 1, a, key);
	return offsearch(mid + 1, top, a, key);
}
int align_paired_read(char * buffer, char *cigar1[20], char *cigar2[20], bwa_seq_t *seqs1, bwa_seq_t *seqs2, hash_element *table1, hash_element *table2,uint64_t level, aryana_args *options,int **d, char **arr, char* tmp_cigar)
{
	int flag1 = 4, flag2=1;
	int ind1=0,ind2 = 0;
	char rnext1[100],rnext2[100];
	long long tlen1=0,tlen2 = 0;
	int head=0;
	//continue;
	int best1[100],best2[100];
	int best_size1=options->potents;
	int best_size2=options->potents;
	int ii=0;
	int jj=0;
	for (ii=0; ii<best_size1; ii++)
	{
		best1[ii]=-1;
		cigar1[ii][0]=0;
	}
	for (ii=0; ii<best_size2; ii++)
	{
		best2[ii]=-1;
		cigar2[ii][0]=0;
	}

	options->len1=seqs1->len;
	options->len2=seqs2->len;

	options->pairID=0;
	int best_num1=best_size1;
	int best_num2=best_size2;
	int best_found1=-1;
	int best_found2=-1;
	long long tlen=0;
	tlen1=0;
	tlen2=0;

	aligner(bwt, seqs1->len, seqs1->seq,  level, table1,best1,best_size1, &best_num1, options);	
	aligner(bwt, seqs2->len, seqs2->seq,  level, table2,best2,best_size2, &best_num2, options);	
	int errors1[20];
	int errors2[20];
	int best_errors=-1;
	double prob=0.0;
	double prob1=0.0;
	double prob2=0.0;
	for (ii=best_num1; ii<best_size1; ii++)
	{
		if (best1[ii]!=-1)
		{
			errors1[ii]=create_cigar(&table1[best1[ii]], cigar1[ii], seqs1->len, seqs1->seq,bwt->seq_len, d, arr, tmp_cigar);
			double scaled_err=scmax*errors1[ii]/(seqs1->len);
			prob1+=pow(base,-scaled_err);
			if (best_found1==-1 || errors1[ii] < errors1[best_found1])
				best_found1=ii;
		}
	}
	for (ii=best_num2; ii<best_size2; ii++)
	{
		if (best2[ii]!=-1)
		{
			errors2[ii]=create_cigar(&table2[best2[ii]], cigar2[ii], seqs2->len, seqs2->seq,bwt->seq_len, d, arr, tmp_cigar);
			double scaled_err=scmax*errors2[ii]/(seqs2->len);
			prob2+=pow(base,-scaled_err);
			if (best_found2==-1 || errors2[ii] < errors2[best_found2])
				best_found2=ii;
		}
	}

	for (ii=best_num1; ii<best_size1; ii++)
		for (jj=best_num2; jj<best_size2; jj++)
		{
			bwtint_t index1=table1[best1[ii]].index;
			bwtint_t index2=table2[best2[jj]].index;
			int rev1=0,rev2=0;
			if (index1>bwt->seq_len / 2)
			{
				index1=(bwt->seq_len / 2) - (index1 - (bwt->seq_len / 2))-seqs1->len;
				rev1=1;
			}
			if (index2>bwt->seq_len / 2)
			{
				index2=(bwt->seq_len / 2) - (index2 - (bwt->seq_len / 2))-seqs2->len;
				rev2=1;
			}

			if (index2>index1)
				tlen=index2-index1+seqs2->len; //TODO
			else
				tlen=index1-index2+seqs1->len; //TODO
			if (tlen < options->min_dis || tlen > options->max_dis)
				continue;

			if (strcmp(options->ori,"ff")==0 && (rev1!=rev2))
				continue;
			if (strcmp(options->ori,"fr")==0 && ((rev1==rev2) || (rev1 && index1 < index2) || (rev2 && index2 < index1)))
				continue;
			if (strcmp(options->ori,"rf")==0 && ((rev1==rev2) || (rev2 && index1 < index2) || (rev1 && index2 < index1)))
				continue;
			double scaled_err1=scmax*errors1[ii]/(seqs1->len);
			double scaled_err2=scmax*errors2[ii]/(seqs2->len);
			prob+=pow(base,-(scaled_err1+scaled_err2));
			if (best_errors==-1 || best_errors>=errors1[ii]+errors2[jj])
			{
				best_errors=errors1[ii]+errors2[jj];
				best_found1=ii;
				best_found2=jj;
			}

		}

	int paired_matched=0;
	int mate1_matched=0;
	int mate2_matched=0;
	seqs1->mapQ=0;
	seqs2->mapQ=0;
	if (best_errors!=-1)
	{
		paired_matched=1;

		double scaled_err=scmax*errors1[best_found1]/(seqs1->len);
		scaled_err+=scmax*errors2[best_found2]/(seqs2->len);
		prob=1.0-pow(base,-scaled_err)/prob;
		if (prob<=0.00000001)
			prob+=0.00000001;
		double mapq=-5*log10(prob);

		seqs1->mapQ=mapq;
		seqs2->mapQ=mapq;

	}

	if (paired_matched || (best_found1!=-1 && options->discordant))
	{
		seqs1->index=table1[best1[best_found1]].index;
		mate1_matched=1;
		
		double scaled_err=scmax*errors1[best_found1]/(seqs1->len);
		prob1=1.0-pow(base,-scaled_err)/prob1;
		if (prob1<=0.00000001)
			prob1+=0.00000001;
		double mapq=-5*log10(prob1);

		seqs1->mapQ+=mapq;
		seqs1->mapQ/=2;
	}
	if (paired_matched || (best_found2!=-1 && options->discordant))
	{
		seqs2->index=table2[best2[best_found2]].index;
		mate2_matched=1;

		double scaled_err=scmax*errors2[best_found2]/(seqs2->len);
		prob2=1.0-pow(base,-scaled_err)/prob2;
		if (prob2<=0.00000001)
			prob2+=0.00000001;
		double mapq=-5*log10(prob2);

		seqs2->mapQ+=mapq;
		seqs2->mapQ/=2;
	}
	flag1=1;
	flag2=1;
	flag1 |= 64;
	flag2 |= 128;
	if (paired_matched)
	{
		flag1 |= 2;
		flag2 |= 2;
	}

	if (mate1_matched==0)
	{
		flag1 |= 4;
		flag2 |= 8;
	}
	else if(seqs1->index <= bwt->seq_len / 2)
		flag1 |= 0;
	else
	{
		flag1 |= 16;
		flag2 |= 32;
	}
	if(flag1 & 16)
		seqs1->index = (bwt->seq_len / 2) - (seqs1->index - (bwt->seq_len / 2))-2*seqs1->len;

	if (mate2_matched==0)
	{
		flag2 |= 4;
		flag1 |= 8;
	}
	else if(seqs2->index <= bwt->seq_len / 2)
		flag2 |= 0;
	else
	{
		flag2 |= 16;
		flag1 |= 32;
	}
	if(flag2 & 16)
		seqs2->index = (bwt->seq_len / 2) - (seqs2->index - (bwt->seq_len / 2))-2*seqs2->len;
	ind1 = offsearch(0, offInd - 1, offset, seqs1->index);
	ind2 = offsearch(0, offInd - 1, offset, seqs2->index);
	seqs1->index -= offset[ind1];
	seqs2->index -= offset[ind2];
	if (mate1_matched && mate2_matched)// && !(/*TODO: no discr...*/1 && !paired_matched))
	{

		if (strcmp(name[ind1],name[ind2])==0)
		{
			strcpy(rnext1,"=");
			strcpy(rnext2,"=");
			if (seqs2->index>seqs1->index)
			{
				tlen1=seqs2->index-seqs1->index+seqs2->len;
				tlen2=-tlen1;
			}
			else
			{
				tlen2=seqs1->index-seqs2->index+seqs1->len;
				tlen1=-tlen2;
			}
		}
		else
		{
			strcpy(rnext1,name[ind2]);
			strcpy(rnext2,name[ind1]);
		}
		head=sam_generator(buffer, seqs1->name, flag1, name[ind1], seqs1->index+1, seqs1->mapQ, cigar1[best_found1], rnext1, seqs2->index+1, tlen1, seqs1->seq, seqs1->qual, seqs1->len);


		//SEQS2
		head+=sam_generator(buffer+head, seqs2->name, flag2, name[ind2], seqs2->index+1, seqs1->mapQ, cigar2[best_found2], rnext2, seqs1->index+1, tlen2, seqs2->seq, seqs2->qual, seqs2->len);
	}
	else if (!mate1_matched && !mate2_matched) 
	{
		head=sam_generator(buffer, seqs1->name, flag1, "*", (uint64_t)(0), (uint32_t)(0), "*", "*", (uint64_t)(0), (long long)(0), seqs1->seq, seqs1->qual, seqs1->len);
		head+=sam_generator(buffer+head, seqs2->name, flag2, "*", (uint64_t)(0), (uint32_t)(0), "*", "*", (uint64_t)(0), (long long)(0), seqs2->seq, seqs2->qual, seqs2->len);
	}
	else {
		if(mate1_matched)
		{
			head=sam_generator(buffer, seqs1->name, flag1, name[ind1], seqs1->index+1, seqs1->mapQ, cigar1[best_found1], "*", (uint64_t)(0), (long long)(0), seqs1->seq, seqs1->qual, seqs1->len);
			head+=sam_generator(buffer+head, seqs2->name, flag2, "*", (uint64_t)(0), (uint32_t)(0), "*", name[ind1], seqs1->index, 0, seqs2->seq, seqs2->qual, seqs2->len);
		}
		else 
		{
			head=sam_generator(buffer, seqs1->name, flag1, "*", (uint64_t)(0), (uint32_t)(0), "*", name[ind2], seqs2->index, 0, seqs1->seq, seqs1->qual, seqs1->len);
			head+=sam_generator(buffer+head, seqs2->name, flag2, name[ind2], seqs2->index+1, seqs1->mapQ, cigar2[best_found2], "*", (uint64_t)(0), (long long)(0), seqs2->seq, seqs2->qual, seqs2->len);
		}
	}
	return head;

	/*			if(true || j % 10 == 0 || j == n_seqs - 1){
				pthread_mutex_lock(&output);
				fputs(buffer, stdout);
				lasttmpsize = 0;
				pthread_mutex_unlock(&output);
				}
				*/

	//lasttmpsize += snprintf(out_buf + lasttmpsize, max_sam_line, "%s", buffer);
	//if(lasttmpsize >= output_buffer * max_sam_line)
	//	fprintf(stderr, "There is a problem\n");
	//lasttmpsize = 0;*/
}
void paired_multiAligner(const gap_opt_t *opt, aryana_args *options){

	//fprintf(stderr,"paired!!!!\n");
//	fprintf(stderr,"min :: %llu max :: %llu \n ", options->min_dis,options->max_dis);
//	fprintf(stderr,"orientation: %s\n ", options->ori);
	bwa_seq_t *seqs1, *seqs2;
	int n_seqs = 0;
	int total_seqs = 0;
	bwtint_t j = 0;
	hash_element *table1,*table2;
	table1=(hash_element *)malloc(TABLESIZE*(sizeof (hash_element)));
	table2=(hash_element *)malloc(TABLESIZE*(sizeof (hash_element)));
	reset_hash(table1);
	reset_hash(table2);
	seqs1 = (bwa_seq_t*)calloc(each_read_size, sizeof(bwa_seq_t));
	seqs2 = (bwa_seq_t*)calloc(each_read_size, sizeof(bwa_seq_t));
	char buffer[10*max_sam_line];
	char *cigar1[100];
	for (j=0; j<options->potents; j++)
		cigar1[j]=(char *)malloc(MAX_CIGAR_SIZE*(sizeof (char)));
	char *cigar2[100];
	for (j=0; j<options->potents; j++)
		cigar2[j]=(char *)malloc(MAX_CIGAR_SIZE*(sizeof (char)));
	buffer[0] = '\0';
	int lasttmpsize = 0;

	int **d=(int **)malloc(MAX_READ_SIZE*(sizeof (int *)));
	char **arr=(char **)malloc(MAX_READ_SIZE*(sizeof (char *)));
	char *tmp_cigar=(char *)malloc(MAX_READ_SIZE*(sizeof (char)));
	for (j=0; j<MAX_READ_SIZE; j++)
	{
		d[j]=(int *)malloc(100*(sizeof (int)));
		arr[j]=(char *)malloc(100*(sizeof (char)));
	}

	while(true){
		pthread_mutex_lock(&input);
		if((seqs1 = bwa_read_seq(ks1, each_read_size, &n_seqs, opt->mode, opt->trim_qual)) == 0){
			pthread_mutex_unlock(&input);
			return;
		}
		seqs2 = bwa_read_seq(ks2, each_read_size, &n_seqs, opt->mode, opt->trim_qual);
		pthread_mutex_unlock(&input);
		lasttmpsize = 0;
		for(j=0;j<n_seqs;j++)
		{
			lasttmpsize+=align_paired_read(buffer+lasttmpsize, cigar1, cigar2, &seqs1[j], &seqs2[j], table1, table2, total_seqs + j + 1, options, d, arr, tmp_cigar);
			if(j % 5 == 0 || j == n_seqs - 1){
				pthread_mutex_lock(&output);
				fputs(buffer, stdout);
				lasttmpsize = 0;
				pthread_mutex_unlock(&output);
			}
		}
		total_seqs += n_seqs;

		bwa_free_read_seq(n_seqs, seqs1);
	}

}

uint64_t reverse_cigar(char *cigar)
{
	int i=0;
	int cigar_size=0;
	for (i=0; cigar[i]; i++)
		cigar_size++;
	char *cigar_tmp=(char *)malloc((cigar_size+50)*(sizeof (char)));
	int num=0;
	char type=cigar[cigar_size-1];
	int offset=1;
	int lasttmpsize=0;
	int length=0;
	for (i=cigar_size-2; i>=0; i--)
	{
		if (cigar[i]>='0' && cigar[i]<='9' )
		{
			num+=offset*(cigar[i]-'0');
			offset*=10;
		}
		else
		{
			lasttmpsize += snprintf(cigar_tmp + lasttmpsize, 10, "%d%c", num,type);
			if (type=='m' || type=='d')
				length+=num;
			type=cigar[i];
			num=0;
			offset=1;
		}
	}
	lasttmpsize += snprintf(cigar_tmp + lasttmpsize, 10, "%d%c", num,type);
	if (type=='m')
		length+=num;
	for (i=0; i<cigar_size; i++)
		cigar[i]=cigar_tmp[i];
	free(cigar_tmp);
	return length;
}

		



	

//////////////////////////////////

int align_read(char * buffer, char *cigar[20], bwa_seq_t *seq, hash_element *table, uint64_t level, aryana_args *options,int **d, char **arr, char* tmp_cigar){
	buffer[0]='\0';
	int flag = 4;
	int ind = 0;
	bwtint_t pnext = 0;
	long long tlen = 0;
	int best_size=options->potents;
	char ** ttt = &buffer;
	char ** ttt2 = ttt;
	int best[100];
	int ii=0;
	for (ii=0; ii<best_size; ii++)
	{
		best[ii]=-1;
		cigar[ii][0]=0;
	}

	//aligner
	int best_num=best_size;
	aligner(bwt, seq->len, seq->seq,  level, table, best,best_size,&best_num, options);	
	int best_errors=-1;
	int best_found=-1;
	double prob=0.0;
	for (ii=best_num; ii<best_size; ii++)
		if (best[ii]!=-1)
		{
	//		fprintf(stderr,"%llu %d\n",table[best[ii]].index,best[ii]);
			int tmp_errors=create_cigar(&table[best[ii]], cigar[ii], seq->len, seq->seq,bwt->seq_len,d,arr,tmp_cigar);

			seq->index=table[best[ii]].index;
			//seq->index = (bwt->seq_len / 2) - (seq->index - (bwt->seq_len / 2))-seq->len;
//			fprintf(stderr,"%d\n",tmp_errors);
			double scaled_err=scmax*tmp_errors/(seq->len);
			prob+=pow(base,-scaled_err);
//			fprintf(stderr,"%d %lf %llu\n",tmp_errors,scaled_err,seq->index);
			if (best_errors==-1 || best_errors>=tmp_errors)
			{

				best_errors=tmp_errors;
				best_found=ii;
			}
		}

//	prob+=pow(5.0,-(10-best_errors/3));
	double scaled_err=scmax*best_errors/(seq->len);
	prob=1.0-pow(base,-scaled_err)/prob;
	if (prob<=0.00000001)
		prob+=0.00000001;
	double mapq=-5*log10(prob);
//	mapq/=1.5;
	seq->mapQ=mapq;
	//	fprintf(stderr,"%lf\n",mapq);
	seq->index=0;

	if (best_found!=-1)
		seq->index=table[best[best_found]].index;

	if (best_found==-1 || seq->index >= bwt->seq_len)
		flag = 4;
	else if(seq->index <= bwt->seq_len / 2)
		flag = 0;
	else
		flag=16;
	if(flag == 16)
	{
		uint64_t ref_aligned_len=reverse_cigar(cigar[best_found]);
		seq->index = (bwt->seq_len / 2) - (seq->index - (bwt->seq_len / 2))-2*ref_aligned_len;
	}
	ind = offsearch(0, offInd - 1, offset, seq->index);

	seq->index -= offset[ind];

	pnext = 0;
	tlen = 0;

	buffer[0]='\0';
//	if (seq->mapQ!=3)
//		return 0;
	if (flag==4)
		return sam_generator(buffer, seq->name, flag, "*", (uint64_t)(0), (uint32_t)(0), "*", "*", (uint64_t)(0), (long long)(0), seq->seq, seq->qual, seq->len);
	else
		return sam_generator(buffer, seq->name, flag, name[ind], seq->index+1, seq->mapQ, cigar[best_found], "*", (uint64_t)(0), (long long)(0), seq->seq, seq->qual, seq->len);

}

void multiAligner(int tid, const gap_opt_t *opt, aryana_args *options){
	fprintf(stderr,"Thread #%d started\n",tid);
	bwa_seq_t *seqs;
	int n_seqs = 0;
	int total_seqs = 0;
	bwtint_t j = 0;
	//fprintf(stderr,"salam\n");
	//fprintf(stderr,"allocate hash\n");
	//hash_element table[4000];
	//	hash_element table[TABLESIZE];
	hash_element *table;
	table=(hash_element *)malloc(TABLESIZE*(sizeof (hash_element)));
	reset_hash(table);
	int **d=(int **)malloc(MAX_READ_SIZE*(sizeof (int *)));
	char **arr=(char **)malloc(MAX_READ_SIZE*(sizeof (char *)));
	char *tmp_cigar=(char *)malloc(MAX_READ_SIZE*(sizeof (char)));
	for (j=0; j<MAX_READ_SIZE; j++){
		d[j]=(int *)malloc(300*(sizeof (int)));
		arr[j]=(char *)malloc(300*(sizeof (char)));
/*		int del=0; 
		for (del=0; del<100; del++)
		{
			d[j][del]=(int *)malloc(5*(sizeof (int)));
			arr[j][del]=(char *)malloc(5*(sizeof (char)));
		}
*/	}
	//fprintf(stderr,"reset hashed\n");
	//fprintf(stderr,"khodafez\n");
	seqs = (bwa_seq_t*)calloc(each_read_size, sizeof(bwa_seq_t));
	//fprintf(stderr,"seqs created\n");
	char *buffer=(char *)malloc(100*max_sam_line*(sizeof (char)));
	buffer[0] = '\0';
	int lasttmpsize = 0;
	char *cigar[100];
	for (j=0; j<100; j++)
		cigar[j]=(char *)malloc(MAX_CIGAR_SIZE*(sizeof (char)));
	//fprintf(stderr, "thread %d starting...\n", tid);
	while(true){
		//fprintf(stderr,"read start %d\n",tid);
		pthread_mutex_lock(&input);
		if((seqs = bwa_read_seq(ks, each_read_size, &n_seqs, opt->mode, opt->trim_qual)) == 0){
			//finish = true;
			pthread_mutex_unlock(&input);
			//fprintf(stderr, "thread %d ending...\n", tid);
			return;
		}
		pthread_mutex_unlock(&input);


		lasttmpsize = 0;
		for(j=0;j<n_seqs;j++){
			lasttmpsize+=align_read(buffer+lasttmpsize, cigar, &seqs[j], table, total_seqs + j + 1, options,d,arr,tmp_cigar);
			if(j % 100 == 0 || j == n_seqs - 1){
				pthread_mutex_lock(&output);
				fputs(buffer, stdout);
				pthread_mutex_unlock(&output);
				lasttmpsize = 0;
			}
		}
		total_seqs += n_seqs;

		bwa_free_read_seq(n_seqs, seqs);
	}
	free(seqs);
	for (j=0; j<MAX_READ_SIZE; j++){
		free(d[j]);
		free(arr[j]);
	}
	free(d);
	free(arr);
	free(table);
	free(tmp_cigar);
	free(buffer);
}

typedef struct {
	int tid;
	const gap_opt_t *opt;
	aryana_args *options;
} thread_args;

void *worker2(void *data){
	//	fprintf(stderr,"worker2\n");
	thread_args *d = (thread_args*)data;
	//	fprintf(stderr,"args done\n");
	//(d->options)->paired=0;
	if ((d->options)->paired==0)
		multiAligner(d->tid, d->opt, d->options);
	else
		paired_multiAligner(d->opt, d->options);
	//	fprintf(stderr,"calling done\n");
	return 0;
}

uint64_t  getIndexInt(uint32_t place){
	uint32_t block = place / (sizeof(uint64_t) * 4);
	int offset = place % (sizeof(uint64_t) * 4);
	uint64_t mask = 3;
	mask = mask & (ref[block] >> (2 * offset));
	return mask;
}

char atom2[4]={'A','C','G','T'};

char getIndexChar(uint32_t place){
	//place += 43;
	uint32_t block = place / (sizeof(uint64_t) * 4);
	int offset = place % (sizeof(uint64_t) * 4);
	//fprintf(stderr, "block = %lld\n", block);
	uint64_t mask = 3;
	mask = mask & (reference[block] >> (2 * offset));
	return atom2[mask];
}

int ref_read(char * file_name){
	fprintf(stderr, "inside ref_read with %s\n", file_name);
	struct stat file_info;
	if(stat(file_name , &file_info) == -1){
		fprintf(stderr, "Could not get the information of file %s\nplease make sure the file exists\n", file_name);
		return -1;
	}
	int fd = open(file_name , O_RDONLY);
	//FILE *fd;
	//fd = xopen(file_name, "rb"); 
	if(fd == -1){
		fprintf(stderr, "Could not open the file %s\nplease make sure the file exists\n", file_name);
		return -1;
	}
	off_t file_size_bytes = file_info.st_size;
	reference_size = ceil ( ((double)file_size_bytes) / (double)(sizeof(uint64_t)) );
	fprintf(stderr, "reference_size = %u\n", reference_size);
	reference_reminder = file_size_bytes % sizeof(uint64_t) ;
	//reference = new base64 [ reference_size ];
	reference = (uint64_t *)malloc(reference_size * sizeof(uint64_t));
	memset ( reference , 0 , reference_size * sizeof(uint64_t) );
	size_t read_size2 = 0;//there is a read_size defined above
	size_t signal;
	size_t total_size = (file_size_bytes);
	unsigned char *uc_buffer = (unsigned char *)(reference);
	int counter=0;

	do{
		signal = read ( fd , (void *)uc_buffer , total_size - read_size2 );
		//signal = fread((void *)uc_buffer, )
		if ( signal == -1 )
		{
			fprintf(stderr, "Error: while writing to file\n");
			if ( close(fd) == -1 )
				fprintf(stderr, "Error: while closing file\n");
			return -1;
		}
		counter++;
		read_size2 += signal;
		uc_buffer += signal;
	}
	while ( read_size2 < total_size );
	if ( close(fd) == -1 )
	{
		fprintf(stderr, "Unable to close the file\n");
		return -1;
	}

	return 0;
}

char getNuc(uint64_t place, uint64_t seq_len){
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
}



//void bwa_aln_core2(const char *prefix, const char *fn_fa, const char *fn_fa1, const char *fn_fa2, const gap_opt_t *opt, pair_opt *options)
void bwa_aln_core2(aryana_args *args)
{
	gap_opt_t opt_tmp;
	opt_tmp.mode = 0;
	opt_tmp.trim_qual = 0;
	gap_opt_t *opt = &opt_tmp;
	char *prefix = args->reference;
	char *fn_fa = args->read_file;
	char *fn_fa1 = args->read_file1;
	char *fn_fa2 = args->read_file2;
	if (args->paired==0)
		fprintf(stderr, "bwa_aln_core2 => prefix: %s, fn_fa: %s\n+", prefix, fn_fa);
	else
		fprintf(stderr, "bwa_aln_core2 => prefix: %s, fn_fa1: %s, fn_fa2: %s\n", prefix, fn_fa1,fn_fa2);
	int i, j;

	fprintf(stderr,"before open...\n");
	if (args->paired==0)
		ks = bwa_open_reads(opt->mode, fn_fa);
	else
	{
		ks1 = bwa_open_reads(opt->mode, fn_fa1);
		ks2 = bwa_open_reads(opt->mode, fn_fa2);
	}
	// initialization


	fprintf(stderr,"loading...\n");
	{ // load BWT
		char *str = (char*)calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
		//strcpy(str, prefix); strcat(str, ".pac");
		//strcpy(str, prefix); strcat(str, ".pac");
		//strcpy(str, "human_g1k_v37.bin");
		strcpy(str, prefix); strcat(str, ".bin");
		//FILE *fp;
		//fp = xopen(str, "rb");
		//int pac_size;
		//pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
		//reference = (bwtint_t*)calloc(pac_size, 1);
		//fread(reference, 1, pac_size, fp);
		ref_read(str);
		free(str);
		//fclose(fp);
	}

	fprintf(stderr,"finished loading\n");

	memset(offset, 0, sizeof(offset));
	char *str = (char*)calloc(strlen(prefix) + 10, 1);
	strcpy(str, prefix); strcat(str, ".ann");
	FILE * ann = fopen(str, "r");
	free(str);

	char line[1000];
	if(fgets(line, sizeof line, ann) == NULL){
		fprintf(stderr, "Error: Empty file\n");
		bwt_destroy(bwt);
		if (args->paired==0)
			bwa_seq_close(ks);
		else
		{
			bwa_seq_close(ks1);
			bwa_seq_close(ks2);
		}
		return;
	}
	fscanf(ann, "%d %s", &i, name[0]);
	while(fgets(line, sizeof line, ann) != NULL){
		fscanf(ann, "%lu", &offset[offInd++]);
		if(fgets(line, sizeof line, ann) == NULL){
			fprintf(stderr, "Error: Empty file\n");
			bwt_destroy(bwt);
			if (args->paired==0)
				bwa_seq_close(ks);
			else
			{
				bwa_seq_close(ks1);
				bwa_seq_close(ks2);
			}
			return;
		}
		fscanf(ann, "%d %s", &i, name[offInd]);
	}

	offset[offInd] = bwt->seq_len / 2;
	char buffer[200000];
	buffer[0] = '\0';
	sam_headers(buffer,offset, offInd);
	//printf("%s", buffer);
	//write(1, buffer, strlen(buffer));
	fputs(buffer, stdout);

	//FILE * fout = fopen("errors.txt", "w");


	fprintf(stderr, "primary:%" PRIu64 "\n", bwt->primary);
	fprintf(stderr, "seq_len: %"PRIu64"\n", bwt->seq_len);
	fprintf(stderr, "bwt_size: %"PRIu64"\n", bwt->bwt_size);
	fprintf(stderr, "sa_intv: %d\n", bwt->sa_intv);
	fprintf(stderr, "n_sa: %"PRIu64"\n", bwt->n_sa);

	pthread_t *threads;
	pthread_attr_t attr;
	thread_args *data;
	threads = (pthread_t*)calloc(args->threads, sizeof(pthread_t));
	pthread_attr_init(&attr);
	data = (thread_args*)calloc(args->threads, sizeof(thread_args));
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	for(j=0;j<args->threads;j++){
		data[j].tid = j;
		data[j].opt = opt;
		data[j].options = args;
	}
	/*	if (options->paired==0)
		singleAligner(opt,options);
		else
		paired_singlealigner(opt,options);
		*/	// core loop
	//err_fwrite(opt, sizeof(gap_opt_t), 1, stdout);
	//----------------------------------------------------------------------Multi thread
	for (j = 0; j < args->threads; j++)
		pthread_create(&threads[j], NULL, worker2, data + j);

	for (j = 0; j < args->threads; ++j)
		pthread_join(threads[j], 0);

	// destroy
	free(data);
	free(threads);
	bwt_destroy(bwt);
	if (args->paired==0)
		bwa_seq_close(ks);
	else
	{
		bwa_seq_close(ks1);
		bwa_seq_close(ks2);
	}

	//fclose(fout);
}
