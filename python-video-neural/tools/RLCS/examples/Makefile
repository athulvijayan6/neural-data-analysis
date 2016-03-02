RLCSINC = ../include
RLCSLIB = ../lib

SHELL  =  /bin/bash
CC     =  gcc
CFLAGS  =   -Wall

all: ParametersRLCSPTR_test SpotVectorPTR_test

ParametersRLCSPTR_test: ParametersRLCSPTR_test.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o ParametersRLCSPTR_test ParametersRLCSPTR_test.c  -lrlcs -lm

SpotVectorPTR_test: SpotVectorPTR_test.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o SpotVectorPTR_test SpotVectorPTR_test.c  -lrlcs -lm


clean:
	rm -rf  ParametersRLCSPTR_test SpotVectorPTR_test
