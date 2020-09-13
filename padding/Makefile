# configurable variables
SIZE=512

# compiler & flags
CC=gcc
CFLAGS=-DSIZE=$(SIZE) -Wall -Wextra -std=c18 -march=native -mtune=native -O3 -I.
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

all: nopad padk padj

nopad: $(BinDir)/nopad.bin
padj: $(BinDir)/padj.bin
padk: $(BinDir)/padk.bin

$(BinDir)/%.bin: $(ObjDir)/common.o $(ObjDir)/%.o
	${MKDIR} $(dir $@)
	$(call Link, $^, $@)
	$(call colorecho, 2, "'$@' build successfully")

$(ObjDir)/%.o: $(SourceDir)/%.c
	${MKDIR} $(dir $@)
	$(call Compile, $<, $@)


run:
	$(BinDir)/nopad.bin
	$(BinDir)/padk.bin
	$(BinDir)/padj.bin

measure:
	${MKDIR} $(MeasureDir)

perf: measure
	lscpu > "$(PERF_TARGET)"
	$(PERF) $(BinDir)/nopad.bin 2>&1 &>> "$(PERF_TARGET)"
	$(PERF) $(BinDir)/padk.bin 2>&1 &>> "$(PERF_TARGET)"
	$(PERF) $(BinDir)/padj.bin 2>&1 &>> "$(PERF_TARGET)"

mem: measure
	lscpu > "$(MEM_TARGET)"
	$(MEM) $(BinDir)/nopad.bin 2>&1 &>> "$(MEM_TARGET)"
	$(MEM) $(BinDir)/padk.bin 2>&1 &>> "$(MEM_TARGET)"
	$(MEM) $(BinDir)/padj.bin 2>&1 &>> "$(MEM_TARGET)"


benchmark: run

clean:
	rm -f $(ObjDir)/*.o $(BinDir)/*.bin

.PHONY: all benchmark clean measure mem nopad padj padk perf run
