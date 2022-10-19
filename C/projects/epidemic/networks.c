#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

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

static int nearest_nb_node=0;
static int nearest_nb_max_degree=0;
static int* nearest_nb_degree=NULL;
static int* nearest_nb_list=NULL;

static void nearest_nb_initialize() {
    nearest_nb_node=100;
    nearest_nb_max_degree=1;

    size_t size1 = nearest_nb_node*nearest_nb_max_degree;
    size_t size2 = nearest_nb_node;
    nearest_nb_list   = (int*)malloc(sizeof(int)*size1);
    nearest_nb_degree = (int*)malloc(sizeof(int)*size2);

    for(int i=0;i<size1;i++) nearest_nb_list[i]=-1;
    for(int i=0;i<size2;i++) nearest_nb_degree[i]=0;
}

static void nearest_nb_realloc(int nnode, int degree) {
    int nb_node=nearest_nb_node;
    int nb_degree=nearest_nb_max_degree;

    int check=0;
    if(nnode>nb_node) {
        nearest_nb_node*=2;
        check=1;
    }
    if(degree>nb_degree) {
        nearest_nb_max_degree=degree;
        check=1;
    }

    if(check) {
        size_t size1 = nearest_nb_node*nearest_nb_max_degree;
        size_t size2 = nearest_nb_node;
        int* temp1 = (int*)malloc(sizeof(int)*size1);
        int* temp2 = (int*)malloc(sizeof(int)*size2);

        for(int i=0;i<size1;i++) temp1[i]=-1;
        for(int i=0;i<size2;i++) temp2[i]=0;
        for(int i_node=0;i_node<nb_node;i_node++) {
            temp2[i_node]=nearest_nb_degree[i_node];
            for(int i_nb=0;i_nb<nb_degree;i_nb++) {
                int i_p = i_node*nb_degree+i_nb;
                int i_n = i_node*nearest_nb_max_degree+i_nb;
                temp1[i_n] = nearest_nb_list[i_p];
            }
        }

        if(nearest_nb_degree!=NULL) free(nearest_nb_degree);
        if(nearest_nb_list!=NULL) free(nearest_nb_list);
        nearest_nb_list = temp1;
        nearest_nb_degree = temp2;
    }
}

static void nearest_nb_append(int i, int j) {
    int i_degree,j_degree;
    int u=i;
    if(j>i) u=j;

    if(i<nearest_nb_node) { 
        i_degree = nearest_nb_degree[i]+1;
    } else {i_degree=1;}

    if(j<nearest_nb_node) {
        j_degree = nearest_nb_degree[j]+1;
    } else {j_degree=1;}

    int u_degree=i_degree;
    if(j_degree>i_degree) u_degree=j_degree;

    nearest_nb_realloc(u,u_degree);

    nearest_nb_list[i*nearest_nb_max_degree+i_degree-1] = j;
    nearest_nb_list[j*nearest_nb_max_degree+j_degree-1] = i;

    nearest_nb_degree[i] = i_degree;
    nearest_nb_degree[j] = j_degree;
}

void nearest_nb_show() {
    for(int i=0;i<nearest_nb_node;i++) {
        if(nearest_nb_degree[i]>0) {
            printf("%d  | ",i);
            for(int j=0;j<nearest_nb_degree[i];j++) {
                printf("%d ",nearest_nb_list[i*nearest_nb_max_degree+j]);
            }
            printf("\n");
        }
    }
}

int nearest_nb_random_assign(int i, gsl_rng* rng) {
    int k = nearest_nb_degree[i];
    double dis=gsl_rng_uniform_pos(rng)*k;

    return nearest_nb_list[i*nearest_nb_max_degree+(int)dis];
}

int nearest_nb_arg_max_degree() {
    int j=0;
    int max_degree=0;

    for(int i=0;i<nearest_nb_node;i++) {
        if(max_degree<nearest_nb_degree[i]) {
            j=i;
            max_degree=nearest_nb_degree[i];
        }
    }

    return j;
}

int* read_edgelist(char* filename, int* nnode, int* nedge) {
    FILE* fp = fopen(filename,"r");
    if(fp==NULL) {
        printf("Can not find the file: %s\n",filename);
        exit(1);
    }

    nearest_nb_initialize();

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
        nearest_nb_append(i,j);
    }

    *nnode = nnode_temp+1;
    fclose(fp);

    //nearest_nb_show();

    return edges;
}

