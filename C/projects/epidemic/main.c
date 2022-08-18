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

int main() {
    char filename[128] = "/home/alan/Works/path_sampling/networks/jupyters/test.edgelist";
    double alpha=0.7;
    unsigned long int seed=39479832;

    int nnode;
    int nedge;
    int* edges = read_edgelist(filename,&nnode,&nedge);


    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    model* m = sis_model_uniform_infection(alpha,nnode,nedge,edges);


    free_model(m);
}
