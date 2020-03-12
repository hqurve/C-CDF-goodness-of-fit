/*
C file to test CDF functions
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <windows.h>
#include "resourcestack.h"
#include "func.h"


struct histogram{
    decimal underflow;
    decimal overflow;
    decimal* values;
    int size;
};

FILE* openFile(char* directory, char*fileName, char* mode){
    char* full = malloc(sizeof(char) * FILENAME_MAX);
    if (full == NULL) return NULL;
    strcpy(full, directory);
    strcat(full, "\\");
    strcat(full, fileName);

    FILE* file = fopen(full, mode);
    free(full);
    return file;
}
void closeFile(void* file){
    fclose(file);
}


int performTests(FILE* output, char* linePrefix, int count, decimal* observedArr, decimal* expectedArr);
int main(int argc, char** argv){
    srand(time(NULL));
    INTIALIZE_RESOURCE_MANAGER;
    if (argc <= 1){
        printf(
            "Parameters:"
            " [output file name]"
            " [sample size]"
            " [min]"
            " [max]"
            " [grain size]"
            PARAM_HINTS
        );
        R_EXIT_SUCCESS;
    }
    if (argc-1 < 5+PARAM_COUNT){
        printf("Missing parameters");
        R_EXIT_FAILURE;
    }
    char* directoryName = argv[1];
    long sampleSize = atoll(argv[2]);
    double min = atof(argv[3]);
    double max = atof(argv[4]);
    double grainSize = atof(argv[5]);
    struct params* params = getParams(argc-5, argv+6);
    
    PUSH_RESOURCE_MESSAGE(params, free, "Params not generated");

    if (sampleSize < 0){
        printf("Sample size is negative");
        R_EXIT_FAILURE
    }else if(min > max){
        printf("Min is more than max");
        R_EXIT_FAILURE
    }else if (grainSize <= 0){
        printf("Invalid grain size");
        R_EXIT_FAILURE
    }
    CreateDirectory(directoryName, NULL);
    
    int grainCount = (max-min)/grainSize;
    grainSize = (max-min)/grainCount;
    //Array for cdf values
    //arr[i] = P(X <= min + grainSize * i)
    decimal* expectedCDFValues = calloc(grainCount+1, sizeof(decimal));
    decimal* sampleCDFValues = calloc(grainCount+1, sizeof(decimal));

    //Array for histogram (includes underflow and overflow)
    //arr[i] = P(min + grainSize * i < X <= min + grainSize * (i+1))
    //underflow = P(X <= min)
    //overflow = P(X > max)
    struct histogram expectedHistogram={
        .underflow = 0,
        .overflow =0,
        .size = grainCount,
        .values = calloc(grainCount, sizeof(decimal))
    };
    struct histogram sampleHistogram={
        .underflow = 0,
        .overflow =0,
        .size = grainCount,
        .values = calloc(grainCount, sizeof(decimal))
    };


    PUSH_RESOURCE_MESSAGE(expectedCDFValues, free, "Unable to allocate expectedCDFValues");
    PUSH_RESOURCE_MESSAGE(sampleCDFValues, free, "Unable to allocate sampleCDFValues");

    PUSH_RESOURCE_MESSAGE(expectedHistogram.values, free, "Unable to allocate expectedHistogram");
    PUSH_RESOURCE_MESSAGE(sampleHistogram.values, free, "Unable to allocate sampleHistogram");



    for (int i = 0; i<grainCount+1; i++){
        expectedCDFValues[i] = CDF(params, min + grainSize * i);
    }

    expectedHistogram.underflow = expectedCDFValues[0];
    for (int i = 0; i< grainCount; i++){
        expectedHistogram.values[i] = (expectedCDFValues[i+1] - expectedCDFValues[i]) / grainSize;
    }
    expectedHistogram.overflow = 1-expectedCDFValues[grainCount]; 


    for (int i =0 ; i<sampleSize; i++){
        decimal value = generateSample(params);
        int pos = ceil((value - min)/grainSize);
        if (pos <= 0) sampleHistogram.underflow++;
        else if(pos <=grainCount) sampleHistogram.values[pos-1]++;
        else sampleHistogram.overflow++;
    }

    sampleHistogram.underflow /= sampleSize;
    sampleHistogram.overflow /= sampleSize;
    for (int i = 0; i < grainCount; i++){
        sampleHistogram.values[i]/=sampleSize;
    }

    sampleCDFValues[0] = sampleHistogram.underflow;
    for (int i = 0; i< grainCount; i++){
        sampleCDFValues[i+1] = sampleCDFValues[i] + sampleHistogram.values[i];
        sampleHistogram.values[i] /= grainSize;
    }



    //FILE *file_samples = openFile(directoryName, "samples.txt", "w");
    FILE *file_summary = openFile(directoryName, "summary.txt", "w");

    FILE *file_histogram = openFile(directoryName, "histogram.txt", "w");

    FILE *file_CDFValues = openFile(directoryName, "CDFValues.txt", "w");
    

    PUSH_RESOURCE_MESSAGE(file_summary, closeFile, "can't open summary.txt");

    PUSH_RESOURCE_MESSAGE(file_histogram, closeFile, "can't open histogram.txt");

    PUSH_RESOURCE_MESSAGE(file_CDFValues, closeFile, "can't open CDFValues.txt");
    

    fprintf(file_histogram, "Expected\tObserved\n");
    fprintf(file_histogram, "%f\t%f\n", (double) expectedHistogram.underflow, (double) sampleHistogram.underflow);
    for(int i =0; i< grainCount; i++){
        fprintf(file_histogram, "%f\t%f\n", (double) expectedHistogram.values[i], (double)sampleHistogram.values[i]);
    }
    fprintf(file_histogram, "%f\t%f\n", (double) expectedHistogram.overflow, (double) sampleHistogram.overflow);

    fprintf(file_CDFValues, "Expected\tObserved\n");
    for(int i=0; i<grainCount+1; i++){
        fprintf(file_CDFValues, "%f\t%f\n", (double) expectedCDFValues[i], (double) sampleCDFValues[i]);
    }
    FREE_RESOURCE(file_histogram);
    FREE_RESOURCE(file_CDFValues);
    


    fprintf(file_summary, "Parameters: ");
    printParams(file_summary, params);
    fprintf(file_summary, "\n");
    fprintf(file_summary, "Number of samples: %d\n", sampleSize);
    
    fprintf(file_summary, "\n");
    fprintf(file_summary, "min: %f\n", (double) min);
    fprintf(file_summary, "max: %f\n", (double)max);
    fprintf(file_summary, "grain size: %f\n", (double) grainSize);
    fprintf(file_summary, "grain count: %d\n", grainCount);
    fprintf(file_summary, "\n");
    
    fprintf(file_summary, "CDF\n");
    performTests(file_summary, "\t", grainCount+1, sampleCDFValues, expectedCDFValues);

    fprintf(file_summary, "\n");

    fprintf(file_summary, "Histogram\n");
    performTests(file_summary, "\t", grainCount, sampleHistogram.values, expectedHistogram.values);

    R_EXIT_SUCCESS
}

int performTests(FILE* output, char* linePrefix, int count, decimal* observedArr, decimal* expectedArr){
    if (count <=0) return 1;
    #define PRINT_LINE(...) fprintf(output, linePrefix); fprintf(output, __VA_ARGS__); fprintf(output, "\n");
    decimal observedTotal = 0, expectedTotal = 0;
    decimal observedSumOfSquares = 0, expectedSumOfSquares = 0;
    decimal residueMin = 0, residueMax = 0;
    decimal residueAverage = 0, residueAbsAverage = 0, residueSumOfSquares = 0;

    decimal standardError = 0;

    decimal rSquaredObserved = 0, rSquaredExpected = 0;

    for (int i =0 ; i<count; i++){
        decimal observed = observedArr[i];
        decimal expected = expectedArr[i];
        decimal residue = expected - observed;

        observedTotal += observed;
        expectedTotal += expected;

        observedSumOfSquares += observed * observed;
        expectedSumOfSquares += expected * expected;

        residueAverage += residue;
        residueAbsAverage += residue <0 ? -residue: residue;
        residueSumOfSquares += residue * residue;

        if (residue > residueMax){
            residueMax = residue;
        }
        if (residue < residueMin){
            residueMin = residue;
        }
    }

    observedSumOfSquares -= observedTotal * observedTotal / count;
    expectedSumOfSquares -= expectedTotal * expectedTotal / count;

    rSquaredObserved = 1 - residueSumOfSquares / observedSumOfSquares;
    rSquaredExpected = 1 - residueSumOfSquares / expectedSumOfSquares;

    residueAverage /= count;
    residueAbsAverage /= count;

    standardError = sqrt(residueSumOfSquares/count);

    PRINT_LINE("count: %d", count);

    PRINT_LINE("Totals:");
    PRINT_LINE("\tExpected: %f", (double) expectedTotal);
    PRINT_LINE("\tObserved: %f", (double) observedTotal);

    PRINT_LINE("Sum of Squares:");
    PRINT_LINE("\tExpected: %f", (double)expectedSumOfSquares);
    PRINT_LINE("\tObserved: %f", (double)observedSumOfSquares);

    PRINT_LINE("Residues:");
    PRINT_LINE("\tMin:    %f", (double) residueMin);
    PRINT_LINE("\tMax:     %f", (double) residueMax);
    PRINT_LINE("\tMean:    %f", (double) residueAverage);
    PRINT_LINE("\t|Mean|:  %f", (double) residueAbsAverage);
    PRINT_LINE("\tSum of Squares: %f", (double) residueSumOfSquares);
    PRINT_LINE("\tStd Error: %f", (double)standardError);

    PRINT_LINE("R-squared:");
    PRINT_LINE("\tExpected: %f", (double) rSquaredExpected);
    PRINT_LINE("\tObserved: %f", (double) rSquaredObserved);

    #undef PRINT_LINE
}