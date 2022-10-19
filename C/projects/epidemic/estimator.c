#include <stdio.h>
#include <stdlib.h>


static int filled_flag=0;
static int sequence_tag=0;
static int sequence_nobs=0;
static int sequence_size=0;
static double* sequence=NULL;
static double* autocorrelation=NULL;

void sequence_malloc(int size, int nobs) {
    sequence = (double*)malloc(sizeof(double)*size*nobs);
    autocorrelation = (double*)malloc(sizeof(double)*size*nobs);

    for(int i=0;i<size;i++){
        for(int j=0;j<nobs;j++) {
            sequence[i*nobs+j]=0;
            autocorrelation[i*nobs+j]=0;
        }
    }

    sequence_size = size;
    sequence_nobs = nobs;
}

static int sequence_append_count=0;
void sequence_append(double* samples) {
    int size = sequence_size;
    int nobs = sequence_nobs;
    int tag  = sequence_tag;
    int index=tag*nobs;

    for(int i=0;i<nobs;i++) {
        sequence[index+i] = samples[i];
    }

    if(filled_flag) {
        for(int i=1;i<(size+1);i++) {
            index = ((tag+i)%size)*nobs;
            for(int j=0;j<nobs;j++) {
                autocorrelation[(size-i)*nobs+j] += samples[j]*sequence[index+j];
            }
        }
        sequence_append_count++;
    }
    sequence_tag++;

    if(sequence_tag==size) {
        filled_flag=1;

        sequence_tag=0;
    }

    if(sequence_append_count==size) {
        FILE* file_a = fopen("autocorrelation.txt","a");
        for(int i=0;i<size;i++) {
            for(int j=0;j<nobs;j++) {
                autocorrelation[i*nobs+j] = autocorrelation[i*nobs+j]/sequence_append_count;
                fprintf(file_a,"%.12e ",autocorrelation[i*nobs+j]);

                autocorrelation[i*nobs+j]=0;
            }
            fprintf(file_a,"\n");
        }
        fclose(file_a);

        sequence_append_count=0;
    }
}
