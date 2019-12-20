/*
C file to test CDF functions
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "resourcestack.h"
#include "func.h"

FILE* openFile(char* directory, char*fileName, char* mode){
    char* full = malloc(sizeof(char) * (strlen(directory) + strlen(fileName)+1));
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

double calcChiSquared(int count, double* observedArray, double* expectedArray){
    double chiSquared = 0;
    for (int i =0; i<count; i++){
        double observed = observedArray[i];
        double expected = expectedArray[i];
        double value = observed - expected;
        value = value * value;
        if (expected != 0)chiSquared += value/ expected;
        else chiSquared += value;
    }
    return chiSquared;
}
int main(int argc, char** argv){
    INTIALIZE_RESOURCE_MANAGER;
    if (argc <= 1){
        printf(
            "Parameters:"
            " [sample size]"
            " [min]"
            " [max]"
            " [grain size]"
            PARAM_HINTS
        );
        R_EXIT_SUCCESS;
    }
    if (argc-1 < 4+PARAM_COUNT){
        printf("Missing parameters");
        R_EXIT_FAILURE;
    }
    char* directoryName = argv[1];
    int sampleSize = atoi(argv[2]);
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

    
    int grainCount = (max-min)/grainSize;
    //Array for cdf values
    //arr[i] = P(X < min + grainSize * i)
    double* expectedCDFValues = malloc(sizeof(double) * grainCount+1);
    double* sampleCDFValues = malloc(sizeof(double) * grainCount+1);
    double* residualCDFValues = malloc(sizeof(double) * grainCount+1);

    //Array for histogram (includes underflow and overflow)
    //arr[0] = P(X < min)
    //arr[i] = P(min + grainSize * (i-1) < X < min + grainSize * i) for i between 1 and grainCount (inclusive)
    //arr[grainCount+1] = P(X > min + grainSize * grainCount)
    double* expectedHistogram = malloc(sizeof(double) * grainCount + 2);
    double* sampleHistogram = malloc(sizeof(double) * grainCount + 2);
    double* residualHistogram = malloc(sizeof(double) * grainCount + 2);


    PUSH_RESOURCE_MESSAGE(expectedCDFValues, free, "Unable to allocate expectedCDFValues");
    PUSH_RESOURCE_MESSAGE(sampleCDFValues, free, "Unable to allocate sampleCDFValues");
    PUSH_RESOURCE_MESSAGE(residualCDFValues, free, "Unable to allocate residualCDFValues");

    PUSH_RESOURCE_MESSAGE(expectedHistogram, free, "Unable to allocate expectedHistogram");
    PUSH_RESOURCE_MESSAGE(sampleHistogram, free, "Unable to allocate sampleHistogram");
    PUSH_RESOURCE_MESSAGE(residualHistogram, free, "Unable to allocate residualHistogram");



    for (int i = 0; i<grainCount+1; i++){
        expectedCDFValues[i] = CDF(params, min + grainSize * i);
    }

    expectedHistogram[0] = sampleSize * expectedCDFValues[0];
    for (int i = 1; i< grainCount+1; i++){
        expectedHistogram[i] = sampleSize * (expectedCDFValues[i] - expectedCDFValues[i-1]);
    }
    expectedHistogram[grainCount+1] = sampleSize*(1-expectedCDFValues[grainCount]); 



    for (int i = 0; i< grainCount+2; i++){
        sampleHistogram[i] = 0;
    }
    for (int i =0 ; i<sampleSize; i++){
        double value = generateSample(params);
        int pos = (value - min)/grainSize;
        if (pos < 0) pos = 0;
        else if(pos <grainCount) pos = pos + 1;
        else pos = grainCount + 1;
        sampleHistogram[pos]++;
    }

    sampleCDFValues[0] = sampleHistogram[0]/sampleSize;
    for (int i = 1; i< grainCount+1; i++){
        sampleCDFValues[i] = sampleCDFValues[i-1] + sampleHistogram[i]/sampleSize;
    }


    //Performing goodness of fit tests
    for (int i =0 ; i<grainCount+1; i++){
        residualCDFValues[i] = sampleCDFValues[i] - expectedCDFValues[i];
    }
    for (int i =0 ; i<grainCount+2; i++){
        residualHistogram[i] = sampleHistogram[i] - expectedHistogram[i];
    }

    double chiSquaredCDF = sampleSize * calcChiSquared(grainCount+1, sampleCDFValues, expectedCDFValues);
    double chiSquaredHistogram = calcChiSquared(grainCount+2, sampleHistogram, expectedHistogram);

    double residualSquaredCDFValues = 0;
    double residualSquaredHistogram = 0;
    for (int i = 0; i<grainCount+1; i++){
        residualSquaredCDFValues += residualCDFValues[i]*residualCDFValues[i];
    }
    residualSquaredCDFValues *= sampleSize;
    residualSquaredCDFValues *= sampleSize;
    for (int i = 0; i<grainCount+2; i++){
        residualSquaredHistogram += residualHistogram[i]*residualHistogram[i];
    }

    //FILE *file_samples = openFile(directoryName, "samples.txt", "w");
    FILE *file_summary = openFile(directoryName, "summary.txt", "w");

    FILE *file_expectedHistogram = openFile(directoryName, "expectedHistogram.txt", "w");
    FILE *file_sampleHistogram = openFile(directoryName, "sampleHistogram.txt", "w");
    FILE *file_residualHistogram = openFile(directoryName, "residualHistogram.txt", "w");

    FILE *file_expectedCDF = openFile(directoryName, "expectedCDF.txt", "w");
    FILE *file_sampleCDF = openFile(directoryName, "sampleCDF.txt", "w");
    FILE *file_residualCDF = openFile(directoryName, "residualCDF.txt", "w");
    

    PUSH_RESOURCE_MESSAGE(file_summary, closeFile, "can't open summary.txt");

    PUSH_RESOURCE_MESSAGE(file_expectedHistogram, closeFile, "can't open expectedHistogram.txt");
    PUSH_RESOURCE_MESSAGE(file_sampleHistogram, closeFile, "can't open sampleHistogram.txt");
    PUSH_RESOURCE_MESSAGE(file_residualHistogram, closeFile, "can't open residualHistogram.txt");

    PUSH_RESOURCE_MESSAGE(file_expectedCDF, closeFile, "can't open expectedCDF.txt");
    PUSH_RESOURCE_MESSAGE(file_sampleCDF, closeFile, "can't open sampleCDF.txt");
    PUSH_RESOURCE_MESSAGE(file_residualCDF, closeFile, "can't open residualCDF.txt");
    

    for(int i =0; i< grainCount+2; i++){
        fprintf(file_expectedHistogram, "%f\n", expectedHistogram[i]);
        fprintf(file_sampleHistogram, "%f\n", sampleHistogram[i]);
        fprintf(file_residualHistogram, "%f\n", residualHistogram[i]);
    }
    for(int i=0; i<grainCount+1; i++){
        fprintf(file_expectedCDF, "%f\n", expectedCDFValues[i]);
        fprintf(file_sampleCDF, "%f\n", sampleCDFValues[i]);
        fprintf(file_residualCDF, "%f\n", residualCDFValues[i]);
    }

    


    fprintf(file_summary, "Parameters: ");
    printParams(file_summary, params);
    fprintf(file_summary, "\n");
    fprintf(file_summary, "Number of samples: %d\n", sampleSize);
    
    fprintf(file_summary, "\n");
    fprintf(file_summary, "min: %f\n", min);
    fprintf(file_summary, "max: %f\n", max);
    fprintf(file_summary, "grain size: %f\n", grainSize);
    fprintf(file_summary, "grain count: %d\n", grainCount);
    fprintf(file_summary, "\n");
    
    fprintf(file_summary, "CDF\n");
    fprintf(file_summary, "\tChi squared: %f\n", chiSquaredCDF);
    fprintf(file_summary, "\tResidual squared: %f\n", residualSquaredCDFValues);
    fprintf(file_summary, "\tResidual standard deviation: %f\n", sqrt(residualSquaredCDFValues/sampleSize));

    fprintf(file_summary, "\n");

    fprintf(file_summary, "Histogram\n");
    fprintf(file_summary, "\tChi squared: %f\n", chiSquaredHistogram);
    fprintf(file_summary, "\tResidual squared: %f\n", residualSquaredHistogram);
    fprintf(file_summary, "\tResidual standard deviation: %f\n", sqrt(residualSquaredHistogram/sampleSize));

    R_EXIT_SUCCESS
}