#include <stdio.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "models.h"
#include "update.h"
#include "stats.h"

int main() {
    int x = 32;
    int y = 32;
    double beta = 10;
    double q3 = 1.5;
    unsigned long int seed = 33474234;
    int thermal = 1000;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    model* m = jq3_ladder_square(x,y,q3);
    world_line* w = malloc_world_line(1024,2*(m->mhnspin),m->nsite);

    estimator* samples[6];
    samples[0] = malloc_estimator(100000,"ms1");
    samples[1] = malloc_estimator(100000,"ms2");
    samples[2] = malloc_estimator(100000,"ms4");
    samples[3] = malloc_estimator(100000,"Xs");
    samples[4] = malloc_estimator(100000,"Xu");
    samples[5] = malloc_estimator(100000,"Noo");

    w->beta = beta;
    for(int i=0;i<(w->nsite);i++) {
        w->istate[i] = 1;
        if(gsl_rng_uniform_pos(rng)<0.5) {
            w->istate[i] = -1;
        }
    }

    for(int i=0;i<thermal;i++) {
        remove_vertices(w);
        insert_vertices(w,m,rng);
        clustering(w,m);
        flip_cluster(w,rng);
    }

    free_model(m);
    free_world_line(w);
    free_estimator(samples[0]);
    free_estimator(samples[1]);
    free_estimator(samples[2]);
    free_estimator(samples[3]);
    free_estimator(samples[4]);
    free_estimator(samples[5]);
    gsl_rng_free(rng);
}
