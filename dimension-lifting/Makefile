# configurable variables
THREADS=4
SIZE=512
NTILES=4
KSPLIT=4

# compiler & flags
CC=gcc
CFLAGS=-DSIZE=$(SIZE) -DTHREADS=$(THREADS) -DNTILES=$(NTILES) \
       -DKSPLIT=$(KSPLIT) -Wall -Wextra -std=c18 -march=native \
       -mtune=native -fopenmp -O3 -I.
LDFLAGS=$(CFLAGS) -lm

MKDIR = -@mkdir -p

# output directories
ObjDir=obj
BinDir=bin
MeasureDir=analysis

# input directories
SourceDir=src

# analysis utils
PERF=perf stat -B -e task-clock,cache-references,cache-misses,cycles,instructions,branches,branch-misses
MEM=/usr/bin/time -v

TARGET := $(shell lscpu | grep 'Model name:' | tr -s ' ' | cut -d':' -f2 | sed 's/^ //' | sed 's/ /_/g')

PERF_TARGET := $(MeasureDir)/$(TARGET).perf
MEM_TARGET := $(MeasureDir)/$(TARGET).mem

# utils
define colorecho
@tput setaf $1
@echo $2
@tput sgr0
endef

define Compile
$(call colorecho, 2, "Compiling object '$@'")
$(CC) -c $< $(CFLAGS) -o $@
endef

define Link
$(call colorecho, 2, "Linking executable '$@'")
$(CC) $1 $(LDFLAGS) -o $@
endef

all: S MDL MDLSL MDLTM

S: $(BinDir)/S.bin
MDL: $(BinDir)/MDL.bin
MDLSL: $(BinDir)/MDLSL.bin
MDLTM: $(BinDir)/MDLTM.bin

$(BinDir)/%.bin: $(ObjDir)/common.o $(ObjDir)/%.o
	${MKDIR} $(dir $@)
	$(call Link, $^, $@)
	$(call colorecho, 2, "'$@' build successfully")

$(ObjDir)/%.o: $(SourceDir)/%.c
	${MKDIR} $(dir $@)
	$(call Compile, $<, $@)


run:
	$(BinDir)/S.bin
	$(BinDir)/MDL.bin
	$(BinDir)/MDLSL.bin
	$(BinDir)/MDLTM.bin

measure:
	${MKDIR} $(MeasureDir)

perf: measure
	lscpu > "$(PERF_TARGET)"
	$(PERF) $(BinDir)/S.bin 2>&1 &>> "$(PERF_TARGET)"
	$(PERF) $(BinDir)/MDL.bin 2>&1 &>> "$(PERF_TARGET)"
	$(PERF) $(BinDir)/MDLSL.bin 2>&1 &>> "$(PERF_TARGET)"
	$(PERF) $(BinDir)/MDLTM.bin 2>&1 &>> "$(PERF_TARGET)"

mem: measure
	lscpu > "$(MEM_TARGET)"
	$(MEM) $(BinDir)/S.bin 2>&1 &>> "$(MEM_TARGET)"
	$(MEM) $(BinDir)/MDL.bin 2>&1 &>> "$(MEM_TARGET)"
	$(MEM) $(BinDir)/MDLSL.bin 2>&1 &>> "$(MEM_TARGET)"
	$(MEM) $(BinDir)/MDLTM.bin 2>&1 &>> "$(MEM_TARGET)"


benchmark: run

clean:
	rm -f $(ObjDir)/*.o $(BinDir)/*.bin

.PHONY: all benchmark clean measure mem S MDL MDLSL MDLTM perf run
