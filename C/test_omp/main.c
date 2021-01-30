#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <omp.h>

double d_t (struct timespec t1, struct timespec t2){
    return (t2.tv_sec-t1.tv_sec)+(double)(t2.tv_nsec-t1.tv_nsec)/1000000000.0;
}

int main() {
    int nt=4;
    int n=100000000;
    int seed = 238723;

    struct timespec t1, t2, t3;

    double* rn1 = (double*)malloc(sizeof(double)*nt*n);
    double* rn2 = (double*)malloc(sizeof(double)*nt*n);

    //serial
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    for(int thread_id=0;thread_id<nt;thread_id++) {
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng,seed*(thread_id+1));

        for(int i=0;i<n;i++) {
            rn1[thread_id*n+i] = sqrt(gsl_rng_uniform_pos(rng));
        }

        gsl_rng_free(rng);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);

    //parallelized
    omp_set_num_threads(nt);
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();    

        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng,seed*(thread_id+1));

        for(int i=0;i<n;i++) {
            rn2[thread_id*n+i] = sqrt(gsl_rng_uniform_pos(rng));
        }

        gsl_rng_free(rng);

    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t3);

    for(int i=0;i<nt*n;i++) {
        double test = rn1[i]-rn2[i];
        if(test!=0){
            printf("%f %f %f\n",test,rn1[i],rn2[i]);
        }
    }

    printf("sequencial for loop:\t %f seconds\n", d_t(t1,t2));
    printf("paralel    for loop:\t %f seconds\n", d_t(t2,t3));
    printf("effection          :\t %f \n",(double)d_t(t1,t2)/(double)d_t(t2,t3));
}
