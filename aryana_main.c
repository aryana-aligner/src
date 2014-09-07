#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "aryana_args.h"
#include "utils.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.6.0-r85"
#endif

int get_value_string(int argc, const char **argv, char * arg, char * value){
	int i = 1;
	for(i=1;i<argc - 1;i++)
		if(strcmp(argv[i], arg) == 0){
			strcpy(value, argv[i + 1]);
			//fprintf(stderr, "value = %s, arg = %s\n", value, arg);
			return 1;
		}
	return 0;
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: bwa (alignment via Burrows-Wheeler transformation)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   bwa <command> [options]\n\n");
	fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
	fprintf(stderr, "         aln           gapped/ungapped alignment\n");
	fprintf(stderr, "         aln2           gapped/ungapped alignment with new algorithm\n");
	fprintf(stderr, "         samse         generate alignment (single ended)\n");
	fprintf(stderr, "         sampe         generate alignment (paired ended)\n");
	fprintf(stderr, "         bwasw         BWA-SW for long queries\n");
	fprintf(stderr, "         fastmap       identify super-maximal exact matches\n");
	fprintf(stderr, "         fa2bin       construcs binary compressed file of fasta\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         pac2bwtgen    alternative algorithm for generating BWT\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n");
	fprintf(stderr, "         pac2cspac     convert PAC to color-space PAC\n");
	fprintf(stderr, "         stdsw         standard SW/NW alignment\n");
	fprintf(stderr, "\n");
	return 1;
}

void bwa_print_sam_PG()
{
	printf("@PG\tID:bwa\tPN:bwa\tVN:%s\n", PACKAGE_VERSION);
}

int main(int argc, char *argv[])
{
	aryana_args args;
	args.discordant=1;
	args.threads=1;
	args.potents=10;
	args.seed_length = -1;
	args.best_factor = -1;
	if(argc < 3){
		fprintf(stderr, "Need more inputs\n");
		return -1;
	}
	fprintf(stderr,"main:\n");
	/*char *reference = (char*)calloc(1000, 1);
	if(!get_value_string(argc, argv, "-x", reference)){
		fprintf(stderr, "Reference genome not specified.\n");
		return -2;
	}
	char *fastq = (char*)calloc(1000, 1);
	if(!get_value_string(argc, argv, "-U", fastq)){
		fprintf(stderr, "Read sequencse file not specified.\n");
		return -2;
	}*/
	//fprintf(stderr, "fastq = %s\n", fastq);
	static struct option long_options[] =
	{
		{"qseq", no_argument, 0, 'Q'},//Reads are QSEQ files.
		{"skip", no_argument, 0, 's'},
		{"qupto", no_argument, 0, 'u'},
		{"trim5", no_argument, 0, '5'},
		{"trim3", no_argument, 0, '3'},
		{"phred33", no_argument, 0, 'M'},
		{"phred64", no_argument, 0, '6'},
		{"very-fast", no_argument, 0, 'V'},
		{"fast", no_argument, 0, 'F'},
		{"sensitive", no_argument, 0, 'Y'},
		{"very-sensitive", no_argument, 0, 'Z'},
		{"dpad", no_argument, 0, 'D'},
		{"gbar", no_argument, 0, 'W'},
		{"ignore-quals", no_argument, 0, 'T'},
		{"nofw", no_argument, 0, 'P'},
		{"norc", no_argument, 0, 'R'},
		{"minins", required_argument, 0, 'I'},
		{"maxins", required_argument, 0, 'X'},
		{"fr", no_argument, 0, '7'},
		{"rf", no_argument, 0, '8'},
		{"ff", no_argument, 0, '9'},
		{"no-mixed", no_argument, 0, 'G'},
		{"no-discordant", no_argument, 0, 0},
		{"dovetail", no_argument, 0, 1},
		{"no-contain", no_argument, 0, 2},
		{"no-overlap", no_argument, 0, 3},
		{"time", no_argument, 0, 't'},
		{"quiet", no_argument, 0, 4},
		{"threads", required_argument, 0, 'p'},
		{"reorder", no_argument, 0, 5},
		{"mm", no_argument, 0, 6},
		{"version", no_argument, 0, 8},
		{"help", no_argument, 0, 9},
		{"seed", required_argument, 0, 10},
		{"factor", required_argument, 0, 'F'}
	//	{0, 0, 0, 0}
	};
	int option_index = 0;
	int c;
	while((c = getopt_long(argc, argv, "x:1:2:U:S:qfrcs:u:5:3:N:L:k:I:X:tp:hP:R:", long_options, &option_index)) >= 0){
		switch(c){
			case 0:
				args.discordant=0;
				break;
			case 1:
				break;
			case 'x':
				args.reference = (char *)malloc(strlen(optarg));
				strcpy(args.reference, optarg);
				break;
			case 'U':
				args.read_file = (char *)malloc(strlen(optarg));
				strcpy(args.read_file, optarg);
				args.single = 1;
				args.paired = 0;
				break;
			case '1':
				args.paired = 1;
				args.read_file1 = (char *)malloc(strlen(optarg));
				strcpy(args.read_file1, optarg);
				break;
			case '2':
				args.paired = 1;
				args.read_file2 = (char *)malloc(strlen(optarg));
				strcpy(args.read_file2, optarg);
				break;
			case '7':
				strcpy(args.ori, "fr");
				break;
			case '8':
				strcpy(args.ori, "rf");
				break;
			case '9':
				strcpy(args.ori, "ff");
				break;
			case 'I':
				args.min_dis = atoi(optarg);
				break;
			case 'X':
				args.max_dis = atoi(optarg);
				break;
			case 'p':
				args.threads = atoi(optarg);
				break;
			case 'P':
				args.potents = atoi(optarg);
				break;
			case 10:
				args.seed_length = atoi(optarg);
				break;
			case 'F':
				args.best_factor = atoi(optarg);
				break;
		}
	}
	bwa_aln_core2(&args);
	//fprintf(stderr, "ori = %s\n", args.ori);
	//bwa_aln_single(args.reference, args.fq);
/*	pair_opt options;
	options.paired=0;
	options.min_dis=0;
	options.max_dis=0;
	if (options.paired==0)
		bwa_aln_core2(argv[optind], argv[optind+1],NULL,NULL, opt, &options);
	else
		bwa_aln_core2(argv[optind], NULL, argv[optind+1], argv[optind+2], opt, &options);*/
	return 0;
}
