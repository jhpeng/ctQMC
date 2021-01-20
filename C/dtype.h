#ifndef dtype_h
#define dtype_h

typedef int (*insert_rule)(int*);

typedef struct model {
    int* bond2type;
    int* bond2hNspin;
    double* bond2weight;
    int** bond2index;
    int** link;
    insert_rule* insert;
    int nsite;
    int nbond;
    double sweight;
} model;

typedef struct vertex {
    double tau;
    int bond;
    int hnspin;
    int* state;
} vertex;

typedef struct world_line {
    vertex* sequenceA;
    vertex* sequenceB;
    int* cluster;
    int* weight;
    int nvertices;
    int mnspin;
    int flag;
    int* state;
    int* last;
    int* first;
    int nsite;
    double beta;
} world_line;

typedef struct estimator {
    char* name;
    double* sample;
    int n;
    double* block;
    int bsize;
    int nblock;
} estimator;

#endif
