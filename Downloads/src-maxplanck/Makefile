CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2
CFLAGS2=	-Wall -O2
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD #-D_NO_SSE2 #-D_FILE_OFFSET_BITS=64
OBJS=		QSufSort.o bwt_gen.o utils.o bwt.o bwtio.o bwtaln.o bwa2.o bwtgap.o sam.o hash.o smith.o aligner.o fa2bin.o \
			is.o bntseq.o bwtmisc.o bwtindex.o ksw.o stdaln.o simple_dp.o \
			bwaseqio.o bwase.o bwape.o kstring.o cs2nt.o \
			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
			bwtsw2_chain.o bamlite.o fastmap.o bwtsw2_pair.o
PROG=		bwa
INCLUDES=	
LIBS=		-lm -lz -lpthread
SUBDIRS=	. bwt_gen

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG) aryana

bwa:$(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

aryana:$(OBJS) aryana_main.o
		$(CC) $(CFLAGS2) $(DFLAGS) $(OBJS) aryana_main.o -o aryana $(LIBS)

QSufSort.o:QSufSort.h

bwt.o:bwt.h
bwtio.o:bwt.h
bwtaln.o:bwt.h bwtaln.h kseq.h bwa2.h
bwt1away.o:bwt.h bwtaln.h
bwt2fmv.o:bwt.h
bntseq.o:bntseq.h
bwtgap.o:bwtgap.h bwtaln.h bwt.h
fastmap:bwt.h
aligner.o: aligner.h bwt.h hash.h smith.h
fa2bin.o: fa2bin.h
hash.o: hash.h
bwa2.o: bwa2.h sam.h aligner.h bwt.h
sam.o: sam.h bwt.h
smith.o: smith.h bwt.h
bwtsw2_core.o:bwtsw2.h bwt.h bwt_lite.h stdaln.h
bwtsw2_aux.o:bwtsw2.h bwt.h bwt_lite.h stdaln.h
bwtsw2_main.o:bwtsw2.h
aryana_main.o: aryana_main.h aryana_args.h bwt.h bwtaln.h kseq.h bwa2.h

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
gen:
	g++ readgen.cpp -Wall -O2 -o readgen
