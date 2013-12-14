#include <string.h>
#include <stdio.h>
#include "fa2bin.h"

int fa2bin(int argc, char *argv[]){
	fprintf(stderr, "inside fa2bin\n");
	if(argc < 2){
		fprintf(stderr, "No fasta specified.\n");
		return -1;
	}
	char * prefix = argv[1];
	char *str = (char*)calloc(strlen(prefix) + 10, 1);
	strcpy(str, prefix);
	fprintf(stderr, "fasta file is: %s\n", str);
	
	FILE *fin = fopen(str, "r");
	strcat(str, ".bin");
	FILE *fbin = fopen(str, "w");
	fprintf(stderr, "output bin is: %s\n", str);

	char c;
	if(fscanf(fin, "%c", &c) == EOF){
		fprintf(stderr, "empty fasta file\n");
		return -2;
	}
	if(c != '>'){
		fprintf(stderr, "wrong fasta format\n");
		return -3;
	}
	else
		while(c != '\n' && fscanf(fin, "%c", &c) != EOF);
		
	int b = 0;
	int i = 0;
	int j = 0;
	while(fscanf(fin, "%c", &c) != EOF){
		if(c == '>'){
			while(c != '\n' && fscanf(fin, "%c", &c) != EOF);
			continue;
		}
		if(c == '\n' || c == '\r')
			continue;
		//if(j && j % 60 == 0)
		//	printf("\n");
		//printf("%c", c);
		//j++;	
		if(j > 99999999)
			break;
		if(c == 'C' || c == 'c')
			b |= (1 << i);
		if(c == 'G' || c == 'g')
			b |= (2 << i);
		if(c == 'T' || c == 't' || c == 'U' || c == 'u')
			b |= (3 << i);
		
		i += 2;
	    if(i == 8){
			fputc(b, fbin);
			b = i = 0;
		}
	}
	if(i) 
		fputc(b, fbin);
	fclose(fbin);
	fclose(fin);
	return 0;
}
