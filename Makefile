CC=icpc
OPTIMIZATION= -ipo -O3 -fargument-noalias -align -ansi-alias -axCORE-AVX2 -restrict -fno-exceptions -fp-model fast=2 -no-prec-div
VEC_REPORT=1
CFLAGS= -ansi -openmp -std=c++11 -Iincludes -static -static-intel $(OPTIMIZATION) -vec-report=$(VEC_REPORT)
CFLAGS_DEBUG=-W -openmp -g -restrict -Wall -ansi -std=c++11  -pedantic -Iincludes -O0 -vec-report=$(VEC_REPORT) -check=uninit -debug all
CFLAGS_MIC=$(CFLAGS) -mmic
LDFLAGS= -openmp $(OPTIMIZATION) -vec-report=$(VEC_REPORT)
LDFLAGS_DEBUG=-openmp -g -vec-report=$(VEC_REPORT) -check=uninit -debug all
LDFLAGS_MIC=-openmp $(OPTIMIZATION) -vec-report=$(VEC_REPORT) -mmic
EXEC=run
SRC=$(wildcard src/*.cpp)
OBJ=$(SRC:src/%.cpp=obj/%.o)
OBJ_MIC=$(SRC:src/%.cpp=obj/%.omic)
OBJ_DEBUG=$(SRC:src/%.cpp=obj/%.odeb)

TEAM_ID = a792ffc4d23771b1447c133*********

all:$(EXEC)

mic:$(OBJ_MIC)
	$(CC) -o $@ $^ $(LDFLAGS_MIC)

default:all

debug:$(OBJ_DEBUG)
	$(CC) -o $@ $^ $(LDFLAGS_DEBUG)


run:$(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

obj/%.o:src/%.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

obj/%.omic:src/%.cpp
	$(CC) -o $@ -c $< $(CFLAGS_MIC)

obj/%.odeb:src/%.cpp
	$(CC) -o $@ -c $< $(CFLAGS_DEBUG)

clean:
	rm -f obj/*.o* $(EXEC) *.zip debug

zip: clean
ifdef TEAM_ID
	zip $(strip $(TEAM_ID)).zip -9r Makefile src/ obj/ includes/ README
else
	@echo "you need to put your TEAM_ID in the Makefile"
endif

submit: zip
ifdef TEAM_ID
	curl -F "file=@$(strip $(TEAM_ID)).zip" -L http://www.intel-software-academic-program.com/contests/ayc/upload/upload.php
else
	@echo "you need to put your TEAM_ID in the Makefile"
endif


