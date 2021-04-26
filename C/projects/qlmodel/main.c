#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <omp.h>

#include "dtype.h"
#include "models.h"
#include "update_omp.h"
#include "stats.h"
#include "union_find.h"

int Lx;
int Ly;
double Lambda;
double Beta;
unsigned long int Seed;

static int referent_conf(int* state){
    if(state[0]*state[1]==1){
        if(state[2]*state[3]==1){
            if(state[1]*state[2]==-1){
                return 1;
            }
        }
    }
    return 0;
}

int count_referent_conf(world_line_omp* w) {
    int count=0;
    for(int y=0;y<Ly;y++) {
        for(int x=0;x<Lx;x++) {
            int u1 = x+y*Lx;
            int u2 = (x+1)%Lx+y*Lx+Lx*Ly;
            int u3 = x+((y+1)%Ly)*Lx;
            int u4 = x+y*Lx+Lx*Ly;

            int state[4];
            state[0] = w->istate[u1];
            state[1] = w->istate[u2];
            state[2] = w->istate[u3];
            state[3] = w->istate[u4];

            if(referent_conf(state)) count++;
        }
    }

    return count;
}

void measurement(world_line_omp* w, model* m) {
    int winding_x = 0;
    int winding_y = 0;

    for(int i=0;i<Lx*Ly;i++) {
        winding_x += w->istate[i];
        winding_y += w->istate[i+Lx*Ly];
    }

    winding_x = winding_x/Lx;
    winding_y = winding_y/Ly;
    printf("winding number x : %d\n",winding_x);
    printf("winding number y : %d\n",winding_y);

    int n_ref_conf = count_referent_conf(w);
    printf("nref : %d \n",n_ref_conf);
}

int main(int argc, char** argv) {
    Lx = 64;
    Ly = 64;
    Lambda = 0.8;
    Beta = 64.0;
    Seed = 389724;

    int nthread = 4;
    int thermal = 1000;

    // Setting up the random number generator
    gsl_rng* rng[nthread];
    for(int i=0;i<nthread;i++){
        rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng[i],(unsigned long int)(Seed*(i+1)*3.14159265359));
    }

    // Setting up the model and the world_line_omp
    model* m = quantum_link_model_2d_square(Lx,Ly,Lambda);

    int mcap = 2.0*(m->sweight)*Beta/nthread;
    if(mcap<1024) mcap=1024;

    world_line_omp* w = malloc_world_line_omp(mcap,2*(m->mhnspin),(m->nsite),nthread);

    w->beta = Beta;

    // Setting up the initial condition
    for(int y=0;y<Ly;y++) {
        for(int x=0;x<Lx;x++) {
            w->istate[x+y*Lx]         = 1-(x+y)%2*2;
            w->istate[x+y*Lx + Lx*Ly] = (x+y)%2*2-1;
/*
            if(y==Ly/2) {
                w->istate[x+y*Lx] = 1;
            }

            if(x==Lx/2){
                w->istate[x+y*Lx+Lx*Ly] = -1;
            }*/
        }
    }
    for(int j=1;j<nthread;j++) {
        for(int i=0;i<(w->nsite);i++) {
            w->istate[i+(w->nsite)*j] = w->istate[i];
        }
    }

    // Thermalization
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

        measurement(w,m);
    }
}
