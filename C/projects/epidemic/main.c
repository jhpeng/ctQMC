#include <stdio.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "sis_models.h"
#include "update.h"

static int* append_edge(int* edges, int* nedge, int i, int j) {
    int buffer=1024;
    int n = *nedge;
    if(n==0) {
        edges = (int*)malloc(sizeof(int)*buffer*2);
    }

    if((n+1)%buffer==0) {
        int* edges_temp = (int*)malloc(sizeof(int)*((n+1)/buffer+1)*buffer*2);
        for(int k=0;k<2*n;k++) {
            edges_temp[k] = edges[k];
        }
        free(edges);
        edges = edges_temp;
    }

    edges[2*n+0] = i;
    edges[2*n+1] = j;
    n++;
    *nedge = n;

    return edges;
}

int* read_edgelist(char* filename, int* nnode, int* nedge) {
    FILE* fp = fopen(filename,"r");
    if(fp==NULL) {
        printf("Error!\n");
        exit(1);
    }

    int* edges = NULL;

    int i,j;
    char data[128];
    int nnode_temp = 0;
    *nedge = 0;
    while(fscanf(fp,"%d %d %s",&i,&j,data)==3) {
        if(nnode_temp<i) {
            nnode_temp=i;
        } else if(nnode_temp<j) {
            nnode_temp=j;
        }
        
        edges = append_edge(edges,nedge,i,j);
    }

    *nnode = nnode_temp+1;
    fclose(fp);

    return edges;
}

void boundary_condition_frozen_initial_state(world_line* w, model* m) {
    int nnode = m->nsite;
    int nbond = m->nbond;
    int length = nnode+(w->nvertices);
    realloc_world_line(w,length);

    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    int n=0;
    for(int i=0;i<nnode;i++) {
        (sequence2[n]).tau      = 0.0;
        (sequence2[n]).bond     = nbond+i;
        (sequence2[n]).hNspin   = 1;
        (sequence2[n]).state[0] = w->istate[i];
        (sequence2[n]).state[1] = w->istate[i];
        n++;
    }
    for(int i=0;i<(w->nvertices);i++) {
        copy_vertex(&(sequence2[n]),&(sequence1[i]));
        n++;
    }

    w->nvertices = n;
    w->flag = !(w->flag);
}

int main() {
    char filename[128] = "/home/alan/Works/path_sampling/networks/jupyters/test.edgelist";
    double alpha=0.7;
    double T = 10.0;
    unsigned long int seed=39479832;

    int thermal = 100000;
    //int nsweep  = 1000;

    int nnode;
    int nedge;
    int* edges = read_edgelist(filename,&nnode,&nedge);


    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    model* m = sis_model_uniform_infection(alpha,nnode,nedge,edges);
    world_line* w = malloc_world_line(1024,2*(m->mhnspin),m->nsite);
    w->beta = T;

    for(int i=0;i<(w->nsite);i++) {
        w->istate[i] = -1;
        if(gsl_rng_uniform_pos(rng)<0.05) {
            w->istate[i] = 1;
        }
    }

    for(int i=0;i<thermal;i++) {
        remove_vertices(w);
        insert_vertices(w,m,rng);
        boundary_condition_frozen_initial_state(w,m);
        clustering(w,m);
        flip_cluster(w,rng);
    }


    free_world_line(w);
    free_model(m);
    gsl_rng_free(rng);
    free(edges);
}
