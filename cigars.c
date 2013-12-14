//YA HAGH
#include<cstdio>

const int READ_SIZE=100;
int match_start[READ_SIZE]={0,4,8,13};
int match_index[READ_SIZE]={2,6,9,14};
int d[READ_SIZE][10];
char arr[READ_SIZE][10];
char tmp_cigar[READ_SIZE];
bool matched[READ_SIZE]={1,1,1,1};
int parts=4;

char read[READ_SIZE]=   "ACCCTAACCCGGA";
char ref[READ_SIZE]="TTACCCGACCCGGA";

/*
char ref[READ_SIZE]="TT ACCC-GAC- CCGGA";
char read[READ_SIZE]=  "ACCC-TAAC-CCGGA";
*/
int smith_waterman(int part,char *cigar,int head)
{
	const int off=5;
	int i=0,j=0;
	for (i=off/2; i<off; i++)
		d[0][i]=i-off/2;
	 
	fprintf(stderr," : %d\n",match_index[part+1]-match_index[part]);
	for (i=1; i<=match_start[part+1]-match_start[part]; i++)
	{
		for (j=0; j<off; j++)
		{
			int real_off=j-off/2;
			int ref_i=i+real_off;
			if (ref_i < 0 )
				continue;
			if (ref_i==0)
			{
				d[i][j]=i;
				continue;
			}

			d[i][j]=d[i-1][j];
			arr[i][j]='m';
//			fprintf(stderr,"%c\n",ref[match_index[part]+ref_i-1]);
			if (ref[match_index[part]+ref_i-1]!=read[match_start[part]+i-1])
			{
	//			fprintf(stderr,"mismatched!\n");
				d[i][j]++;
			}
			if (j>0 && (d[i][j] > d[i][j-1]+1) )
			{
				d[i][j]=d[i][j-1]+1;
				arr[i][j]='d'; 
			}
			if (j<off-1 && (d[i][j] > d[i-1][j+1]+1))
			{
				d[i][j]=d[i-1][j+1]+1;
				arr[i][j]='i'; 		
			}
			fprintf(stderr,"i=%d, ref_i=%d, d=%d, arr=%c\n",i,ref_i,d[i][j],arr[i][j]);
		}
	}
/* do ta tike darim mikhaym match konim,
   d[i][j]= read[match[start]] ta khode read[match[start]+i] ro ba
   ref[match[index]] ta khode ref[match[index]+ i + j
   */
	int cur_off=(match_index[part+1]-match_index[part])-(match_start[part+1]-match_start[part])+off/2;
	int cur_i=match_start[part+1]-match_start[part];
/*	int best_d=d[match_start[part+1]-match_start[part]][0];

	for (i=1; i<off; i++)
	{
			fprintf(stderr," ::: %d\n",d[match_start[part+1]-match_start[part]][i]);
		if (d[match_start[part+1]-match_start[part]][i]<best_d)
		{
			best_d=d[match_start[part+1]-match_start[part]][i];
			cur_off=i;
		}
	}
*/	fprintf(stderr,"curf off : %d\n",cur_off);

	fprintf(stderr,"d[ans] : %d\n",d[match_start[part+1]-match_start[part]][cur_off]);
	int tail=0;
	while(true)
	{
		int ref_i=cur_i+cur_off-off/2;
		if (cur_i==0)
		{
			for (j=0; j<ref_i; j++)
				tmp_cigar[tail++]='d';
			break;
		}
		if (ref_i==0)
		{

			fprintf(stderr,"heey %d %d\n",cur_i,cur_off);
			for (j=0; j<cur_i; j++)
				tmp_cigar[tail++]='i';
			break;
		}
		tmp_cigar[tail++]=arr[cur_i][cur_off];
		if (arr[cur_i][cur_off]=='i')
		{
			cur_off++;
			cur_i--;
		}	
		else if (arr[cur_i][cur_off]=='d')
		{
			cur_off--;
		}
		else if (arr[cur_i][cur_off]=='m')
		{
			cur_i--;
		}
		fprintf(stderr,"bug %d %d\n",cur_i,cur_off);
		fprintf(stderr,"%c",tmp_cigar[tail-1]);
	}
	tmp_cigar[tail]=0;
	fprintf(stderr,"%s\n",tmp_cigar);
	int ct=0;
	char last=tmp_cigar[tail-1];
	for (i=tail-1; i>=-1; i--)
	{
		
		fprintf(stderr,"-- %d\n",tail);
		
		if (i>=0 && tmp_cigar[i]==last )
			ct++;
		else
		{
			head+=snprintf(cigar+head,4,"%d%c",ct,last);
			last=tmp_cigar[i];
			ct=1;
		}
	}
	printf("%s\n",cigar);
	return head;
}
// -----+++++++
//*----***++++++
/*
char ref[READ_SIZE]="TT ACCC-GAC- CCGGA";
char read[READ_SIZE]=  "ACCC-TAAC-CCGGA";
*/
void build_cigar(char * cigar)
{
	int i=0;
	int head=0;
	fprintf(stderr,"%d\n",parts);
	for (i=0; i<parts-1; i++) //TODO:check if reverse and if reversed check all the code
	{
		if (matched[i])
		{
/*			int shift=0;
			if (i>0 && matched[i-1])
				shift=(match_index[i] - match_index[i-1]) - (match_start[i]-match_start[i-1]);
			fprintf(stderr,"shift = %d %d %d\n",shift,i,match_index[i]);
			if (shift>0)
				head+=snprintf(cigar+head,4,"%di",shift);
			if (shift<0)
				head+=snprintf(cigar+head,4,"%dd",-shift);*/
			int ref_len=match_index[i+1]-match_index[i];
			int read_len=match_start[i+1]-match_start[i];
			int min_len=ref_len;
			if (read_len<min_len)
				min_len=read_len;
			head+=snprintf(cigar+head,4,"%dm", min_len);
			if (ref_len<read_len)
				head+=snprintf(cigar+head,4,"%di", read_len-ref_len);
			if (ref_len>read_len)
				head+=snprintf(cigar+head,4,"%dd", ref_len-read_len);
		}
		else
		{
			head=smith_waterman(i,cigar,head);
		}
	}
}
			
int main()
{
	char cigar[200];
	cigar[0]=0;
	build_cigar(cigar);
	printf("%s\n",cigar);
	return 0;
}
