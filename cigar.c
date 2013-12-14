//YA HAGH
#include<cstdio>

const int READ_SIZE=100;
int match_start[READ_SIZE]={0,4,8,13};
int match_index[READ_SIZE]={2,6,9,14};
int d[READ_SIZE][10];
char arr[READ_SIZE][10];
char tmp_cigar[READ_SIZE];
bool matched[READ_SIZE]={1,0,1,1};
int parts=4;

char ref[READ_SIZE]=   "ACCCTAACCCGGA";
char read[READ_SIZE]="TTACCCGACCCGGA";


int smith_waterman(int part,char *cigar,int head)
{
	const int off=5;
	int i=0,j=0;
	for (i=0; i<off/2; i++)
		d[0][i]=off/2-i;
	for (i=off/2; i<off; i++)
		d[0][i]=i-off/2;
	for (i=1; i<=match_start[part+1]-match_start[part]; i++)
	{
		for (j=0; j<off; j++)
		{
			d[i][j]=d[i-1][j];
			arr[i][j]='m';
			if (ref[match_index[part]+i+j-off/2]!=read[match_start[part]+i])
			{
				d[i][j]++;
			}
			if (j>0 && d[i][j] > d[i-1][j-1]+1)
			{
				d[i][j]=d[i-1][j-1]+1;
				arr[i][j]='i'; 
			}
			if (j<off-1 && d[i][j] > d[i-1][j+1]+1)
			{
				d[i][j]=d[i-1][j+1]+1;
				arr[i][j]='d'; 		
			}
		}
	}
/* do ta tike darim mikhaym match konim,
   d[i][j]= read[match[start]] ta khode read[match[start]+i] ro ba
   ref[match[index]] ta khode ref[match[index]+ i + j
   */
	int cur_off=(match_start[part+1]-match_start[part])-(match_index[part+1]-match_index[part])+off/2;
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
*/	fprintf(stderr,"%d\n",cur_off);
	for (i=match_start[part+1]-match_start[part]; i>0; i--)
	{
		tmp_cigar[i]=arr[i][cur_off];
		if (arr[i][cur_off]=='i')
			cur_off--;
		else if (arr[i][cur_off]=='d')
			cur_off++;
		fprintf(stderr,"%c",tmp_cigar[i]);
	}
	fprintf(stderr,"\n");
	int ct=0;
	char last=tmp_cigar[1];
	for (i=1; i<=match_start[part+1]-match_start[part]+1; i++)
	{
		
		if (i<=match_start[part+1]-match_start[part] && tmp_cigar[i]==last )
			ct++;
		else
		{
			head+=snprintf(cigar+head,4,"%d%c",ct,last);
			last=tmp_cigar[i];
			ct=0;
		}
	}
	return head;
}
// -----+++++++
//*----***++++++
void build_cigar(char * cigar)
{
	int i=0;
	int head=0;
	fprintf(stderr,"%d\n",parts);
	for (i=0; i<parts-1; i++) //TODO:check if reverse and if reversed check all the code
	{
		if (matched[i])
		{
			int shift=0;
			if (i>0)
				shift=(match_index[i] - match_index[i-1]) - (match_start[i]-match_start[i-1]);
			fprintf(stderr,"%d\n",shift);
			if (shift>0)
				head+=snprintf(cigar+head,4,"%di",shift);
			if (shift<0)
				head+=snprintf(cigar+head,4,"%dd",-shift);
			head+=snprintf(cigar+head,4,"%dm", (match_start[i+1] - match_start[i]) + ( shift < 0 ? shift : 0));
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
