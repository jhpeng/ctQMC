#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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

static void print_state(int* state, int nnode) {
    for(int i=0;i<nnode;i++) {
        printf("%d ",(state[i]+1)/2);
    }
    printf("\n");
}

void show_configuration(world_line* w, model* m, double* time_list, int ntime) {
    int* pstate = w->pstate;
    int nnode = w->nsite;
    int mhnspin = m->mhnspin;

    for(int i=0;i<nnode;i++) pstate[i] = w->istate[i];

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    int i=0;
    int index, i_node;
    vertex* v;
    for(int n=0;n<(w->nvertices) && i<ntime;n++) {
        v = &(sequence[n]);
        if(time_list[i]<(v->tau)) {
            print_state(pstate,nnode);
            i++;
        }

        for(i_node=0;i_node<(v->hNspin);i_node++) {
            index = m->bond2index[(v->bond)*mhnspin+i_node];
            pstate[index] = v->state[(v->hNspin)+i_node];
        }
    }
    for(;i<ntime;i++) {
        print_state(pstate,nnode);
    }
}

time_t start_time,end_time;
unsigned long int measurement_count=0;
double* infected_ratio=NULL;
void measurement(world_line* w, model* m, double* time_list, int ntime, int block_size) {
    if(infected_ratio==NULL) {
        infected_ratio = (double*)malloc(sizeof(double)*ntime);
        start_time = clock();
    }
    int* pstate = w->pstate;
    int nnode = w->nsite;
    int mhnspin = m->mhnspin;

    for(int i=0;i<nnode;i++) pstate[i] = w->istate[i];

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    int i=0;
    int index, i_node;
    vertex* v;
    for(int n=0;n<(w->nvertices) && i<ntime;n++) {
        v = &(sequence[n]);
        if(time_list[i]<(v->tau)) {
            double ir=0;
            for(int i_node=0;i_node<nnode;i_node++) {
                ir+=0.5*(pstate[i_node]+1);
            }
            ir = ir/nnode;
            infected_ratio[i] += ir;
            i++;
        }

        for(i_node=0;i_node<(v->hNspin);i_node++) {
            index = m->bond2index[(v->bond)*mhnspin+i_node];
            pstate[index] = v->state[(v->hNspin)+i_node];
        }
    }
    for(;i<ntime;i++) {
        double ir=0;
        for(int i_node=0;i_node<nnode;i_node++) {
            ir+=0.5*(pstate[i_node]+1);
        }
        ir = ir/nnode;
        infected_ratio[i] += ir;
    }
    measurement_count++;

    if(measurement_count==block_size) {
        end_time = clock();
        printf("------------------------------\n");
        printf(" t    |    I/N\n");
        for(i=0;i<ntime;i++) {
            infected_ratio[i] = infected_ratio[i]/block_size;
            printf("%.4lf  %.12lf\n",time_list[i]*w->beta,infected_ratio[i]);

            infected_ratio[i] = 0;
        }
        measurement_count=0;

        printf("time for this block = %.2lf(sec)\n",(double)(end_time-start_time)/CLOCKS_PER_SEC);
        start_time = clock();
    }
}

int main() {
    char filename[128] = "/home/alan/Works/path_sampling/networks/jupyters/test.edgelist";
    double alpha=1.0;
    double T = 40.0;
    unsigned long int seed=79636232;

    int block_size=1000;
    int thermal = 1000;
    //int nsweep  = 1000;

    int nnode;
    int nedge;
    int* edges = read_edgelist(filename,&nnode,&nedge);


    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);


    model* m = sis_model_uniform_infection(alpha,nnode,nedge,edges);
    world_line* w = malloc_world_line(1024,2*(m->mhnspin),m->nsite);
    w->beta = T;

    // initial state
    for(int i=0;i<(w->nsite);i++) {
        w->istate[i] = -1;
        //if(gsl_rng_uniform_pos(rng)<0.05) {
        if(i<10) {
            w->istate[i] = 1;
        }
    }

    // thermalization
    for(int i=0;i<thermal;i++) {
        remove_vertices(w);
        swapping_graphs(w,m,rng);
        insert_vertices(w,m,rng);
        boundary_condition_frozen_initial_state(w,m);
        clustering(w,m);
        flip_cluster(w,rng);
    }

    // measurement
    double dt = 2.0;
    int ntime = (int)(T/dt+1);
    double* time_list = (double*)malloc(sizeof(double)*ntime);
    for(int i=0;i<ntime;i++) {
        time_list[i] = (dt*i)/T;
    }

    for(;;) {
        for(int i=0;i<10;i++) {
            remove_vertices(w);
            swapping_graphs(w,m,rng);
            insert_vertices(w,m,rng);
            boundary_condition_frozen_initial_state(w,m);
            clustering(w,m);
            flip_cluster(w,rng);
        }

        measurement(w,m,time_list,ntime,block_size);
    }

    // free memory
    free(infected_ratio);
    free(time_list);
    free_world_line(w);
    free_model(m);
    gsl_rng_free(rng);
    free(edges);
}
