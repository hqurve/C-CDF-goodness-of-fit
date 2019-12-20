#ifndef __FUNC_H__
#define __FUNC_H__

#include <stdio.h> 
#include <stdlib.h>

struct params{
    int i;
};

#define PARAM_COUNT 0
#define PARAM_HINTS ""

struct params* getParams(int argc, char** argv){
    struct params* params = malloc(sizeof(struct params));
    if (params == NULL) return params;

    return params;
}

void printParams(FILE* file, struct params* params){

}
double generateSample(struct params* params){
    return ((double)rand() / RAND_MAX);
}
double CDF(struct params* params, double x){
    return x;
};

#endif