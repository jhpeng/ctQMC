#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <omp.h>

#include "dtype.h"
#include "models.h"
#include "update_omp.h"
#include "stats.h"

int Mode;
int Lx;
int Ly;
double Q3;
double Beta;
unsigned long int Seed;

double* staggered_structure=NULL;
void measurement(world_line_omp* w, model* m, estimator** samples) {
    if(staggered_structure==NULL) {
        staggered_structure = (double*)malloc(sizeof(double)*(Lx*Ly));
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
        w->pstate[i] = w->istate[i];
    }

    int l=0;
    for(int i_thread=0;i_thread<(w->nthread);i_thread++) {

        vertex* sequence = w->sequenceB[i_thread];
        if(w->flag[i_thread])
            sequence = w->sequenceA[i_thread];

        for(int i=0;i<(w->len[i_thread]);i++) {
            int bond   = (sequence[i]).bond;
            int hNspin = (sequence[i]).hNspin;
        
            for(int j=0;j<hNspin;j++) {
                int i_site = m->bond2index[bond*mhnspin+j];
                int statep = (sequence[i]).state[j];
                int staten = (sequence[i]).state[j+hNspin];
                int dif    = (statep*staten-1)*statep;
                ms += staggered_structure[i_site]*dif;
                w->pstate[i_site] = staten;
            }

            msx += ms;
            ms1 += abs(ms);
            ms2 += ms*ms;
            ms4 += ms*ms*ms*ms;
        }

        l += w->len[i_thread];
    }


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

int main(int argc, char** argv) {
    Mode = atoi(argv[1]);
    Lx = atoi(argv[2]);
    Ly = atoi(argv[3]);
    Q3 = atof(argv[4]);
    Beta = atof(argv[5]);
    Seed = atoi(argv[6]);

    int thermal = atoi(argv[7]);
    int nsweep  = atoi(argv[8]);
    int nthread = atoi(argv[9]);

    gsl_rng* rng[nthread];
    for(int i=0;i<nthread;i++) {
        rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng[i],(unsigned long int)(Seed*(i+1)*3.14159265359));
    }

    model* m = jq3_ladder_square(Lx,Ly,Q3);
    if(Mode==1) {
        free_model(m);
        m = jq3_ladder_square_impurity_spin_half(Lx,Ly,Q3);
    } else if(Mode==2) {
        free_model(m);
        m = jq3_ladder_square_impurity_spin_one(Lx,Ly,Q3);
    }

    int mcap = 2.0*(m->sweight)*Beta/nthread;
    if(mcap<1024) mcap=1024;

    world_line_omp* w = malloc_world_line_omp(mcap,2*(m->mhnspin),(m->nsite),nthread);

    int nobs=7;
    estimator* samples[nobs];
    samples[0] = malloc_estimator(nsweep,"ms1");
    samples[1] = malloc_estimator(nsweep,"ms2");
    samples[2] = malloc_estimator(nsweep,"ms4");
    samples[3] = malloc_estimator(nsweep,"Xs");
    samples[4] = malloc_estimator(nsweep,"Xu");
    samples[5] = malloc_estimator(nsweep,"Noo");
    samples[6] = malloc_estimator(nsweep,"time");

    w->beta = Beta;
    for(int i=0;i<(w->nsite);i++) {
        w->istate[i] = 1;
        if(gsl_rng_uniform_pos(rng[0])<0.5) {
            w->istate[i] = -1;
        }

        for(int j=1;j<nthread;j++) 
            w->istate[i+(w->nsite)*j] = w->istate[i];
    }

    double times[50];
    for(int i=0;i<thermal;i++) {
        double start = omp_get_wtime();
        remove_vertices_omp(w);
        double t1 = omp_get_wtime();
        insert_vertices_omp(w,m,rng);
        double t2 = omp_get_wtime();
        clustering_inner_omp(w,m);
        double t3 = omp_get_wtime();
        clustering_crossing_omp(w);
        double t4 = omp_get_wtime();
        flip_cluster_omp(w,rng);
        double end = omp_get_wtime();
        //check_world_line_omp_configuration(w,m);

        int noo=0;
        for(int j=0;j<nthread;j++)
            noo+=w->len[j];

        times[i%50] = end-start;

        printf("-----------------------------------------\n");
        printf("thremal:%d | Noo:%d | nthread=%d ",i,noo,nthread);
        if(i>50) {
            double mtime=0;
            for(int j=0;j<50;j++) mtime+=times[j];
            mtime = mtime/50;
            printf("| time:%f \n",mtime);
        } else {
            printf("\n");
        }
        printf("remove_vertices_omp     : %f  \n",(t1-start)/(end-start));
        printf("insert_vertices_omp     : %f  \n",(t2-t1)/(end-start));
        printf("clustering_inner_omp    : %f  \n",(t3-t2)/(end-start));
        printf("clustering_crossing_omp : %f  \n",(t4-t3)/(end-start));
        printf("flip_cluster_omp        : %f  \n",(end-t4)/(end-start));
    }

    int i=0;
    int nblock=0;
    double timer_print = omp_get_wtime();
    while(i<nsweep) {
        double start = omp_get_wtime();
        remove_vertices_omp(w);
        insert_vertices_omp(w,m,rng);
        clustering_inner_omp(w,m);
        clustering_crossing_omp(w);
        flip_cluster_omp(w,rng);
        measurement(w,m,samples);
        double end = omp_get_wtime();
        append_estimator(samples[6],end-start);

        i++;
        
        if(i>1024 && nblock!=(samples[0])->nblock) {
            if((omp_get_wtime()-timer_print)>60.0) {
                nblock = (samples[0])->nblock;
                printf("=========================================\n");
                printf("L=%d q=%.4f beta=%.1f\n",Lx,Q3,Beta);
                print_detail(samples[1]);
                print_detail(samples[3]);
                print_detail(samples[4]);
                print_detail(samples[5]);
                print_detail(samples[6]);

                for(int i_obs=0;i_obs<nobs;i_obs++)
                    save_estimator(samples[i_obs]);

                timer_print = omp_get_wtime();
            }
        }
    }

    for(int i_obs=0;i_obs<nobs;i_obs++)
        save_estimator(samples[i_obs]);

    free_model(m);
    free_world_line_omp(w);
    free_estimator(samples[0]);
    free_estimator(samples[1]);
    free_estimator(samples[2]);
    free_estimator(samples[3]);
    free_estimator(samples[4]);
    free_estimator(samples[5]);
    free_estimator(samples[6]);

    for(int i=0;i<nthread;i++)
        gsl_rng_free(rng[i]);
}
