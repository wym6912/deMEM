CC?=gcc
CPP = g++

STATIC_LINK = 1

ifeq ($(STATIC_LINK), 1)
	STATIC_ARG = -static
else
	STATIC_ARG = 
endif
WFLAGS= -Wall -Wextra  -DNDEBUG -Wno-ignored-qualifiers
#WFLAGS=-Wall -Werror -Wno-unused-function
OPT_FLAGS= -g -O3 -ffast-math -funroll-loops -m64 -fomit-frame-pointer -fopenmp -D_FILE_OFFSET_BITS=64
CFLAGS = $(WFLAGS) $(OPT_FLAGS) $(STATIC_ARG)
CPPFLAGS = ${CFLAGS} -std=c++17
SYSTEM?=$(shell uname -s)

PREFIX = /usr/local
BINDIR = $(PREFIX)/bin
HAS_64 = 1

THREADS = 1

#valgrind
#CFLAGS += -g -O0
##

MALLOC_COUNT=
MALLOC_COUNT64=

LIBOBJ_SPAWN = lib/ext_prog_call.o lib/gaps.o lib/irreg_part.o lib/seqio.o  \
               lib/libfile_v2.o lib/program_spawn.o lib/ext_prog_config.o lib/parblockaligner.o

LIBOBJ_MAIN = lib/program_spawn.o lib/ext_prog_config.o lib/libfile_v2.o

LIBSSW = external/ssw/src/ssw.o

SUBWMSA_TARGETS = profilealign fragalign profile_two_align

ABPOA_TRAGETS = abpoa abpoa_profile

SUBMAFFT_DIR = external/submafft/Linux/

ifeq ($(HAS_64), 1)
	TARGETS = block_dispatch samemalign_main seq_remove_char splitmem
else
	TARGETS = block_dispatch samemalign_main seq_remove_char splitmem
endif
##

TARGETS_ALL = $(TARGETS) subwmsa $(ABPOA_TRAGETS)

DEBUG = 0
M64 = 0
DNA = 0
TERMINATOR = 1
DEFINES = -DDEBUG=$(DEBUG)
ifeq ($(DNA),1)
	DEFINES += -DDNA=$(DNA)
endif
DEFINES += -DTERMINATOR=$(TERMINATOR)
CFLAGS += $(DEFINES)
ifeq ($(DEBUG), 1)
	CFLAGS += -g -O0
endif
##

##

all:${TARGETS_ALL}

clean:
	${RM} *.o external/*.o lib/*.o external/malloc_count/*.o ${TARGETS} $(SUBWMSA_TARGETS) $(ABPOA_TRAGETS)
	make -C external/ssw/src/ clean
	make -C external/subwmsa clean
	make -C external/abpoa clean

install: all
	mkdir -p $(BINDIR)
	cp $(TARGETS) $(SUBWMSA_TARGETS) $(BINDIR)
	cp -rf $(SUBMAFFT_DIR)/* $(PREFIX)
	@echo Install to $(BINDIR) done.


block_dispatch: block_dispatch.cpp lib/libblockaligner.o $(LIBOBJ_SPAWN)
	$(CPP) -o $@ block_dispatch.cpp ${LIBOBJ_SPAWN} $(LDFLAGS) $(CPPFLAGS) lib/libblockaligner.o $(LIBSSW) $(STATIC_ARG) -Wl,--whole-archive -ldl -pthread -Wl,--no-whole-archive

samemalign_main: samemalign_main.cpp ${LIBOBJ_MAIN}
	$(CPP) -o $@ samemalign_main.cpp ${LIBOBJ_MAIN} $(LDFLAGS) $(CPPFLAGS) $(STATIC_ARG)

ssw:
	make -C external/ssw/src/ -j$(THREADS) ssw.o

subwmsa:
	make -C external/subwmsa -j$(THREADS) $(SUBWMSA_TARGETS)
	@for filename in $(SUBWMSA_TARGETS); do cp external/subwmsa/$$filename .; done

abpoa:
	make -C external/abpoa -j$(THREADS) $(ABPOA_TRAGETS)
	@for filename in $(ABPOA_TRAGETS); do cp external/abpoa/bin/$$filename .; done

lib/libblockaligner.o: lib/libblockaligner.cpp ssw
	$(CPP) $(CPPFLAGS) -c -o $@ $< $(LIBSSW)

lib/parblockaligner.o: lib/parblockaligner.cpp ssw lib/libblockaligner.O
	$(CPP) $(CPPFLAGS) -c -o $@ $< $(LIBSSW) -pthread

seq_remove_char: seq_remove_char.cpp
	$(CPP) -o $@ -static seq_remove_char.cpp

splitmem:
	make -C external/splitmem/src -j$(THREADS)
	mv external/splitmem/src/SplitMEM splitmem
