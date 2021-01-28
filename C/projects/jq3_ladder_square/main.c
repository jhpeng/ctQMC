#include <stdio.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "models.h"
#include "update.h"
#include "stats.h"

int Lx;
int Ly;
double Q3;
double Beta;
unsigned long int Seed;

double* staggered_structure=NULL;
void measurement(world_line* w, model* m, estimator** samples) {
    if(staggered_structure==NULL) {
        staggered_structure = (double*)malloc(sizeof(double)*w->nsite);
        for(int j=0;j<Ly;j++) {
            for(int i=0;i<Lx;i++) {
                int i_site = i+j*Lx;
                staggered_structure[i_site] = ((i+j)%2)*2-1;
            }
        }
    }

    int mhnspin = m->mhnspin;

    double mz = 0;
    double ms = 0;
    double msx = 0;
    double ms1 = 0;
    double ms2 = 0;
    double ms4 = 0;
    double chiu = 0;

    int nsite = w->nsite;
    for(int i=0;i<nsite;i++) {
        mz += w->istate[i];
        ms += (w->istate[i])*staggered_structure[i];
    }

    vertex* sequence = w->sequenceB;
    if(w->flag)
        sequence = w->sequenceA;

    for(int i=0;i<(w->nvertices);i++) {
        int bond   = (sequence[i]).bond;
        int hNspin = (sequence[i]).hNspin;

        for(int j=0;j<hNspin;j++) {
            int i_site = m->bond2index[bond*mhnspin+j];
            int statep = (sequence[i]).state[j];
            int staten = (sequence[i]).state[j+hNspin];
            int dif    = (statep*staten-1)*statep;
            ms += staggered_structure[i_site]*dif;
        }

        msx += ms;
        ms1 += abs(ms);
        ms2 += ms*ms;
        ms4 += ms*ms*ms*ms;
    }

    int l   = w->nvertices;

    chiu = mz*mz*Beta*0.25/nsite;
    if(l!=0) {
        msx = Beta*(msx*msx+ms2)*0.25/l/(l+1)/nsite;
        ms1 = ms1*0.5/nsite/l;
        ms2 = ms2*0.25/nsite/nsite/l;
        ms4 = ms4*0.0625/nsite/nsite/nsite/nsite/l;
    } else {
        msx = Beta*ms*ms/2.0/nsite*0.25;
        ms1 = abs(ms)*0.5/nsite;
        ms2 = ms*ms*0.25/nsite/nsite;
        ms4 = ms*ms*ms*ms/0.0625/nsite/nsite/nsite/nsite;
    }

    append_estimator(samples[0],ms1);
    append_estimator(samples[1],ms2);
    append_estimator(samples[2],ms4);
    append_estimator(samples[3],msx);
    append_estimator(samples[4],chiu);
    append_estimator(samples[5],l);
}

int main() {
    Lx = 256;
    Ly = 256;
    Beta = 20;
    Q3 = 1.5;
    Seed = 3347423;

    int thermal = 10000;
    int nsweep  = 100000;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,Seed);

    model* m = jq3_ladder_square(Lx,Ly,Q3);
    world_line* w = malloc_world_line(1024,2*(m->mhnspin),m->nsite);

    estimator* samples[6];
    samples[0] = malloc_estimator(100000,"ms1");
    samples[1] = malloc_estimator(100000,"ms2");
    samples[2] = malloc_estimator(100000,"ms4");
    samples[3] = malloc_estimator(100000,"Xs");
    samples[4] = malloc_estimator(100000,"Xu");
    samples[5] = malloc_estimator(100000,"Noo");

    w->beta = Beta;
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

        printf("%d \n",w->nvertices);
        if(!check_periodic(w,m)){
            printf("fail of check periodic!\n");
            exit(-1);
        }
    }

    int i=0;
    while(i<nsweep) {
        remove_vertices(w);
        insert_vertices(w,m,rng);
        clustering(w,m);
        flip_cluster(w,rng);
        measurement(w,m,samples);
        i++;

        if(i%1000==0){
            printf("=========================================\n");
            print_detail(samples[1]);
            print_detail(samples[3]);
            print_detail(samples[4]);
            print_detail(samples[5]);
        }
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
