#ifndef __FUNC_H__
#define __FUNC_H__

#include <stdio.h> 
#include <stdlib.h>

/*
    For illustration a uniform distribution is shown
*/
struct params{
    int n;
};

#define PARAM_COUNT 1
#define PARAM_HINTS "[n]"

struct params* getParams(int argc, char** argv){
    struct params* params = malloc(sizeof(struct params));
    if (params == NULL) return params;
    params->n = atoi(argv[0]);
    if (params->n <1){
        free(params);
        return NULL;
    }
    return params;
}

void printParams(FILE* file, struct params* params){
    fprintf(file, "n=%d", params->n);
}
double generateSample(struct params* params){
    return ((double)rand() / RAND_MAX) * params->n;
}
double CDF(struct params* params, double x){
    x /= params-> n;
    if (x <0) return 0;
    if (x >=1) return 1;
    return x;
};

#endif