#include <stdio.h>
#include <stdlib.h>

#include "dtype.h"

static void create_cmf(double* cmf, double* weight, int length) {
    int i=0;

    cmf[0] = weight[i];
    for(i=1;i<length;i++) {
        cmf[i] = cmf[i-1]+weight[i];
    }
}

/* graph name : singlet_proj_1
//     2    3    
//     |....|
//      ....
//     |    |
//     0    1
*/

int link_rule_singlet_proj_1[] = {0,0,2,2,2,1,2,1};

int insert_rule_singlet_proj_1(int* state) {
    if(state[0]*state[1]==-1) {
        return 1;
    }
    return 0;
}

/* graph name : triplet_proj_1
//     2   3    
//      \ /
//       X
//      / \
//     0   1
*/

int link_rule_triplet_proj_1[] = {0,1,1,0,2,2,1,1};

int insert_rule_triplet_proj_1(int* state) {
    if(state[0]*state[1]==1) {
        return 1;
    }
    return 0;
}

/* graph name : singlet_proj_3
//     6    7 8    9 10   11
//     |....| |....| |....|
//      ....   ....   ....
//     |    | |    | |    |
//     0    1 2    3 4    5
*/

int link_rule_singlet_proj_3[] 
        = {0,0,2,2,4,4,6,6,8,8,10,10,2,1,2,1,2,1,2,1,2,1,2,1};

int insert_rule_singlet_proj_3(int* state) {
    if(state[0]*state[1]==-1) {
        if(state[2]*state[3]==-1) {
            if(state[4]*state[5]==-1) {
                return 1;
            }
        }
    }
    return 0;
}

/* graph name : singlet_triplet_proj_2_1
//     6    7 8    9  10  11
//     |....| |....|   \ /
//      ....   ....     X
//     |    | |    |   / \
//     0    1 2    3  4   5
*/

int link_rule_singlet_triplet_proj_2_1[] 
        = {0,0,2,2,4,5,6,6,8,8,5,4,2,1,2,1,2,2,2,1,2,1,1,1};

int insert_rule_singlet_triplet_proj_2_1(int* state) {
    if(state[0]*state[1]==-1) {
        if(state[2]*state[3]==-1) {
            if(state[4]*state[5]==1) {
                return 1;
            }
        }
    }
    return 0;
}

/* graph name : spatial_2_2
//     4    5  6    7
//     |    |  |    |
//     |----|  |----|
//     |    |  |    |
//     0    1  2    3
*/

int link_rule_spatial_2_2[]
        = {0,0,2,2,0,0,2,2,4,1,4,1,1,1,1,1};

/* graph name : spatial_4
//     4    5    6    7
//     |    |    |    |
//     |----|----|----|
//     |    |    |    |
//     0    1    2    3
*/

int link_rule_spatial_4[]
        = {0,0,0,0,0,0,0,0,8,1,1,1,1,1,1,1};

/* graph name : single_box_vertex_1
//     4    5    6    7
//     |____|....|____|
//      ____ .... ____
//     |    |    |    |
//     0    1    2    3 
*/

int link_rule_single_box_vertex_1[]
        = {0,0,0,0,4,4,4,4,4,1,1,1,4,1,1,1};

int insert_rule_single_box_vertex_1(int* state){
    if(state[0]*state[1]==1){
        if(state[2]*state[3]==1){
            if(state[1]*state[2]==-1){
                return 1;
            }
        }
    }
    return 0;
}

/* graph name : single_box_vertex_2
//     4    5    6    7
//     |____|....|____|
//     |    |    |    |
//     0    1    2    3 
*/

int link_rule_single_box_vertex_2[]
        = {0,0,0,0,0,0,0,0,8,1,1,1,1,1,1,1};

int insert_rule_single_box_vertex_2(int* state){
    if(state[0]*state[1]==1){
        if(state[2]*state[3]==1){
            if(state[1]*state[2]==-1){
                return 1;
            }
        }
    }
    return 0;
}

/* graph name : single_box_non_vertex
//     4    5    6    7
//     |    |    |    |
//     |----|----|----|
//     |    |    |    |
//     0    1    2    3 
*/

int link_rule_single_box_non_vertex[]
        = {0,0,0,0,0,0,0,0,8,1,1,1,1,1,1,1};

int insert_rule_single_box_non_vertex(int* state){
    if(state[0]*state[1]==1){
        if(state[2]*state[3]==1){
            if(state[1]*state[2]==-1){
                return 0;
            }
        }
    }
    return 1;
}

model* jq3_ladder_square(int lx, int ly, double q) {
    int nsite = lx*ly;
    int nbond = 4*nsite;
    int mhnspin = 6;
    model* m = malloc_model(nsite,nbond,mhnspin);

    int n=0;
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = x+y*lx;
            int j = (x+1)%lx+y*lx;

            m->bond2type[n]   = 0;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0.5;
            m->sweight += 0.5;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            m->bond2index[n*mhnspin+2] = -1;
            m->bond2index[n*mhnspin+3] = -1;
            m->bond2index[n*mhnspin+4] = -1;
            m->bond2index[n*mhnspin+5] = -1;
            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = x+y*lx;
            int j = x+((y+1)%ly)*lx;

            m->bond2type[n]   = 0;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0.5;
            m->sweight += 0.5;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            m->bond2index[n*mhnspin+2] = -1;
            m->bond2index[n*mhnspin+3] = -1;
            m->bond2index[n*mhnspin+4] = -1;
            m->bond2index[n*mhnspin+5] = -1;
            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i1 = (x+0)%lx+((y+0)%ly)*lx;
            int i2 = (x+1)%lx+((y+0)%ly)*lx;
            int i3 = (x+0)%lx+((y+1)%ly)*lx;
            int i4 = (x+1)%lx+((y+1)%ly)*lx;
            int i5 = (x+0)%lx+((y+2)%ly)*lx;
            int i6 = (x+1)%lx+((y+2)%ly)*lx;

            m->bond2type[n]   = 1;
            m->bond2hNspin[n] = 6;
            m->bond2weight[n] = 0.125*q;
            m->sweight += 0.125*q;
            m->bond2index[n*mhnspin+0] = i1;
            m->bond2index[n*mhnspin+1] = i2;
            m->bond2index[n*mhnspin+2] = i3;
            m->bond2index[n*mhnspin+3] = i4;
            m->bond2index[n*mhnspin+4] = i5;
            m->bond2index[n*mhnspin+5] = i6;
            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i1 = (x+0)%lx+((y+0)%ly)*lx;
            int i2 = (x+0)%lx+((y+1)%ly)*lx;
            int i3 = (x+1)%lx+((y+0)%ly)*lx;
            int i4 = (x+1)%lx+((y+1)%ly)*lx;
            int i5 = (x+2)%lx+((y+0)%ly)*lx;
            int i6 = (x+2)%lx+((y+1)%ly)*lx;

            m->bond2type[n]   = 1;
            m->bond2hNspin[n] = 6;
            m->bond2weight[n] = 0.125*q;
            m->sweight += 0.125*q;
            m->bond2index[n*mhnspin+0] = i1;
            m->bond2index[n*mhnspin+1] = i2;
            m->bond2index[n*mhnspin+2] = i3;
            m->bond2index[n*mhnspin+3] = i4;
            m->bond2index[n*mhnspin+4] = i5;
            m->bond2index[n*mhnspin+5] = i6;
            n++;
        }
    }

    for(int i=0;i<8;i++) {
        m->link[i] = link_rule_singlet_proj_1[i];
    }
    for(int i=0;i<24;i++) {
        m->link[4*mhnspin+i] = link_rule_singlet_proj_3[i];
    }

    m->insert[0] = insert_rule_singlet_proj_1;
    m->insert[1] = insert_rule_singlet_proj_3;

    m->nsite = nsite;
    m->nbond = nbond;
    m->mhnspin = mhnspin;
    create_cmf(m->cmf,m->bond2weight,nbond);

    return m;
}

model* jq3_ladder_square_impurity_spin_half(int lx, int ly, double q) {
    int nsite = lx*ly;
    int nbond = 4*nsite;
    int mhnspin = 6;
    model* m = malloc_model(nsite,nbond,mhnspin);

    int n=0;
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = x+y*lx;
            int j = (x+1)%lx+y*lx;

            m->bond2type[n]   = 0;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0.5;
            m->sweight += 0.5;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            m->bond2index[n*mhnspin+2] = -1;
            m->bond2index[n*mhnspin+3] = -1;
            m->bond2index[n*mhnspin+4] = -1;
            m->bond2index[n*mhnspin+5] = -1;
            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = x+y*lx;
            int j = x+((y+1)%ly)*lx;

            m->bond2type[n]   = 0;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0.5;
            m->sweight += 0.5;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            m->bond2index[n*mhnspin+2] = -1;
            m->bond2index[n*mhnspin+3] = -1;
            m->bond2index[n*mhnspin+4] = -1;
            m->bond2index[n*mhnspin+5] = -1;
            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i1 = (x+0)%lx+((y+0)%ly)*lx;
            int i2 = (x+1)%lx+((y+0)%ly)*lx;
            int i3 = (x+0)%lx+((y+1)%ly)*lx;
            int i4 = (x+1)%lx+((y+1)%ly)*lx;
            int i5 = (x+0)%lx+((y+2)%ly)*lx;
            int i6 = (x+1)%lx+((y+2)%ly)*lx;

            m->bond2type[n]   = 1;
            m->bond2hNspin[n] = 6;
            m->bond2weight[n] = 0.125*q;
            m->sweight += 0.125*q;
            m->bond2index[n*mhnspin+0] = i1;
            m->bond2index[n*mhnspin+1] = i2;
            m->bond2index[n*mhnspin+2] = i3;
            m->bond2index[n*mhnspin+3] = i4;
            m->bond2index[n*mhnspin+4] = i5;
            m->bond2index[n*mhnspin+5] = i6;
            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i1 = (x+0)%lx+((y+0)%ly)*lx;
            int i2 = (x+0)%lx+((y+1)%ly)*lx;
            int i3 = (x+1)%lx+((y+0)%ly)*lx;
            int i4 = (x+1)%lx+((y+1)%ly)*lx;
            int i5 = (x+2)%lx+((y+0)%ly)*lx;
            int i6 = (x+2)%lx+((y+1)%ly)*lx;

            m->bond2type[n]   = 1;
            m->bond2hNspin[n] = 6;
            m->bond2weight[n] = 0.125*q;
            m->sweight += 0.125*q;
            m->bond2index[n*mhnspin+0] = i1;
            m->bond2index[n*mhnspin+1] = i2;
            m->bond2index[n*mhnspin+2] = i3;
            m->bond2index[n*mhnspin+3] = i4;
            m->bond2index[n*mhnspin+4] = i5;
            m->bond2index[n*mhnspin+5] = i6;
            n++;
        }
    }

    for(int i=0;i<8;i++) {
        m->link[i] = link_rule_singlet_proj_1[i];
    }
    for(int i=0;i<24;i++) {
        m->link[4*mhnspin+i] = link_rule_singlet_proj_3[i];
    }

    m->insert[0] = insert_rule_singlet_proj_1;
    m->insert[1] = insert_rule_singlet_proj_3;

    for(int i=0;i<nbond;i++) {
        for(int j=0;j<mhnspin;j++) {
            if((m->bond2index[i*mhnspin+j])==(nsite-1)){
                m->sweight -= m->bond2weight[i];
                m->bond2weight[i]=0;
            }
        }
    }

    m->nsite = nsite-1;
    m->nbond = nbond;
    m->mhnspin = mhnspin;
    create_cmf(m->cmf,m->bond2weight,nbond);

    return m;
}

model* jq3_ladder_square_impurity_spin_one(int lx, int ly, double q) {
    int nsite = lx*ly;
    int nbond = 4*nsite;
    int mhnspin = 6;
    model* m = malloc_model(nsite,nbond,mhnspin);

    int n=0;
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = x+y*lx;
            int j = (x+1)%lx+y*lx;

            m->bond2type[n]   = 0;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0.5;
            m->sweight += 0.5;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            m->bond2index[n*mhnspin+2] = -1;
            m->bond2index[n*mhnspin+3] = -1;
            m->bond2index[n*mhnspin+4] = -1;
            m->bond2index[n*mhnspin+5] = -1;

            if(i==(nsite-1) || j==(nsite-1)) {
                m->bond2type[n] = 2;
            }

            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = x+y*lx;
            int j = x+((y+1)%ly)*lx;

            m->bond2type[n]   = 0;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0.5;
            m->sweight += 0.5;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            m->bond2index[n*mhnspin+2] = -1;
            m->bond2index[n*mhnspin+3] = -1;
            m->bond2index[n*mhnspin+4] = -1;
            m->bond2index[n*mhnspin+5] = -1;

            if(i==(nsite-1) || j==(nsite-1)) {
                m->bond2type[n] = 2;
            }

            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i1 = (x+0)%lx+((y+0)%ly)*lx;
            int i2 = (x+1)%lx+((y+0)%ly)*lx;
            int i3 = (x+0)%lx+((y+1)%ly)*lx;
            int i4 = (x+1)%lx+((y+1)%ly)*lx;
            int i5 = (x+0)%lx+((y+2)%ly)*lx;
            int i6 = (x+1)%lx+((y+2)%ly)*lx;

            m->bond2type[n]   = 1;
            m->bond2hNspin[n] = 6;
            m->bond2weight[n] = 0.125*q;
            m->sweight += 0.125*q;
            m->bond2index[n*mhnspin+0] = i1;
            m->bond2index[n*mhnspin+1] = i2;
            m->bond2index[n*mhnspin+2] = i3;
            m->bond2index[n*mhnspin+3] = i4;
            m->bond2index[n*mhnspin+4] = i5;
            m->bond2index[n*mhnspin+5] = i6;

            if(i1==(nsite-1) || i2==(nsite-1)) {
                m->bond2type[n] = 3;
                m->bond2index[n*mhnspin+0] = i3;
                m->bond2index[n*mhnspin+1] = i4;
                m->bond2index[n*mhnspin+2] = i5;
                m->bond2index[n*mhnspin+3] = i6;
                m->bond2index[n*mhnspin+4] = i1;
                m->bond2index[n*mhnspin+5] = i2;
            } else if(i3==(nsite-1) || i4==(nsite-1)) {
                m->bond2type[n] = 3;
                m->bond2index[n*mhnspin+0] = i1;
                m->bond2index[n*mhnspin+1] = i2;
                m->bond2index[n*mhnspin+2] = i5;
                m->bond2index[n*mhnspin+3] = i6;
                m->bond2index[n*mhnspin+4] = i3;
                m->bond2index[n*mhnspin+5] = i4;
            } else if(i5==(nsite-1) || i6==(nsite-1)) {
                m->bond2type[n] = 3;
            }

            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i1 = (x+0)%lx+((y+0)%ly)*lx;
            int i2 = (x+0)%lx+((y+1)%ly)*lx;
            int i3 = (x+1)%lx+((y+0)%ly)*lx;
            int i4 = (x+1)%lx+((y+1)%ly)*lx;
            int i5 = (x+2)%lx+((y+0)%ly)*lx;
            int i6 = (x+2)%lx+((y+1)%ly)*lx;

            m->bond2type[n]   = 1;
            m->bond2hNspin[n] = 6;
            m->bond2weight[n] = 0.125*q;
            m->sweight += 0.125*q;
            m->bond2index[n*mhnspin+0] = i1;
            m->bond2index[n*mhnspin+1] = i2;
            m->bond2index[n*mhnspin+2] = i3;
            m->bond2index[n*mhnspin+3] = i4;
            m->bond2index[n*mhnspin+4] = i5;
            m->bond2index[n*mhnspin+5] = i6;

            if(i1==(nsite-1) || i2==(nsite-1)) {
                m->bond2type[n] = 3;
                m->bond2index[n*mhnspin+0] = i3;
                m->bond2index[n*mhnspin+1] = i4;
                m->bond2index[n*mhnspin+2] = i5;
                m->bond2index[n*mhnspin+3] = i6;
                m->bond2index[n*mhnspin+4] = i1;
                m->bond2index[n*mhnspin+5] = i2;
            } else if(i3==(nsite-1) || i4==(nsite-1)) {
                m->bond2type[n] = 3;
                m->bond2index[n*mhnspin+0] = i1;
                m->bond2index[n*mhnspin+1] = i2;
                m->bond2index[n*mhnspin+2] = i5;
                m->bond2index[n*mhnspin+3] = i6;
                m->bond2index[n*mhnspin+4] = i3;
                m->bond2index[n*mhnspin+5] = i4;
            } else if(i5==(nsite-1) || i6==(nsite-1)) {
                m->bond2type[n] = 3;
            }

            n++;
        }
    }

    for(int i=0;i<8;i++) {
        m->link[i] = link_rule_singlet_proj_1[i];
    }
    for(int i=0;i<24;i++) {
        m->link[1*4*mhnspin+i] = link_rule_singlet_proj_3[i];
    }
    for(int i=0;i<8;i++) {
        m->link[2*4*mhnspin+i] = link_rule_triplet_proj_1[i];
    }
    for(int i=0;i<24;i++) {
        m->link[3*4*mhnspin+i] = link_rule_singlet_triplet_proj_2_1[i];
    }

    m->insert[0] = insert_rule_singlet_proj_1;
    m->insert[1] = insert_rule_singlet_proj_3;
    m->insert[2] = insert_rule_triplet_proj_1;
    m->insert[3] = insert_rule_singlet_triplet_proj_2_1;

    m->nsite = nsite;
    m->nbond = nbond;
    m->mhnspin = mhnspin;
    create_cmf(m->cmf,m->bond2weight,nbond);

    return m;
}

model* quantum_link_model_2d_square(int lx, int ly, double lambda) {
    int nsite = 2*lx*ly;
    int nbond = 3*lx*ly;
    int mhnspin = 4;
    model* m = malloc_model(nsite,nbond+4*lx*ly,mhnspin);

    int n=0;
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int u1 = x+y*lx;
            int u2 = (x+1)%lx+y*lx+lx*ly;
            int u3 = x+((y+1)%ly)*lx;
            int u4 = x+y*lx+lx*ly;

            m->bond2type[n]   = 0;
            m->bond2hNspin[n] = 4;
            m->bond2weight[n] = 1.0;
            m->sweight += 1.0;
            m->bond2index[n*mhnspin+0] = u1;
            m->bond2index[n*mhnspin+1] = u2;
            m->bond2index[n*mhnspin+2] = u3;
            m->bond2index[n*mhnspin+3] = u4;
            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int u1 = x+y*lx;
            int u2 = (x+1)%lx+y*lx+lx*ly;
            int u3 = x+((y+1)%ly)*lx;
            int u4 = x+y*lx+lx*ly;

            m->bond2type[n]   = 1;
            m->bond2hNspin[n] = 4;
            m->bond2weight[n] = lambda;
            m->sweight += lambda;
            m->bond2index[n*mhnspin+0] = u1;
            m->bond2index[n*mhnspin+1] = u2;
            m->bond2index[n*mhnspin+2] = u3;
            m->bond2index[n*mhnspin+3] = u4;
            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int u1 = x+y*lx;
            int u2 = (x+1)%lx+y*lx+lx*ly;
            int u3 = x+((y+1)%ly)*lx;
            int u4 = x+y*lx+lx*ly;

            m->bond2type[n]   = 2;
            m->bond2hNspin[n] = 4;
            m->bond2weight[n] = 1.0;
            m->sweight += 1.0;
            m->bond2index[n*mhnspin+0] = u1;
            m->bond2index[n*mhnspin+1] = u2;
            m->bond2index[n*mhnspin+2] = u3;
            m->bond2index[n*mhnspin+3] = u4;
            n++;
        }
    }
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int u1 = x+y*lx;
            int u2 = x+y*lx+lx*ly;
            int u3 = (x+lx-1)%lx+y*lx;
            int u4 = x+((y+ly-1)%ly)*lx+lx+ly;

            m->bond2type[n]   = 3;
            m->bond2hNspin[n] = 4;
            m->bond2weight[n] = 0.0;
            m->sweight += 0.0;
            m->bond2index[n*mhnspin+0] = u1;
            m->bond2index[n*mhnspin+1] = u4;
            m->bond2index[n*mhnspin+2] = u2;
            m->bond2index[n*mhnspin+3] = u3;
            n++;

            m->bond2type[n]   = 3;
            m->bond2hNspin[n] = 4;
            m->bond2weight[n] = 0.0;
            m->sweight += 0.0;
            m->bond2index[n*mhnspin+0] = u1;
            m->bond2index[n*mhnspin+1] = u2;
            m->bond2index[n*mhnspin+2] = u3;
            m->bond2index[n*mhnspin+3] = u4;
            n++;

            m->bond2type[n]   = 3;
            m->bond2hNspin[n] = 4;
            m->bond2weight[n] = 0.0;
            m->sweight += 0.0;
            m->bond2index[n*mhnspin+0] = u1;
            m->bond2index[n*mhnspin+1] = u3;
            m->bond2index[n*mhnspin+2] = u2;
            m->bond2index[n*mhnspin+3] = u4;
            n++;

            m->bond2type[n]   = 4;
            m->bond2hNspin[n] = 4;
            m->bond2weight[n] = 0.0;
            m->sweight += 0.0;
            m->bond2index[n*mhnspin+0] = u1;
            m->bond2index[n*mhnspin+1] = u2;
            m->bond2index[n*mhnspin+2] = u3;
            m->bond2index[n*mhnspin+3] = u4;
            n++;
        }
    }

    for(int i=0;i<4*mhnspin;i++) {
        m->link[0*4*mhnspin+i] = link_rule_single_box_vertex_1[i];
    }
    for(int i=0;i<4*mhnspin;i++) {
        m->link[1*4*mhnspin+i] = link_rule_single_box_vertex_2[i];
    }
    for(int i=0;i<4*mhnspin;i++) {
        m->link[2*4*mhnspin+i] = link_rule_single_box_non_vertex[i];
    }
    for(int i=0;i<4*mhnspin;i++) {
        m->link[3*4*mhnspin+i] = link_rule_spatial_2_2[i];
    }
    for(int i=0;i<4*mhnspin;i++) {
        m->link[4*4*mhnspin+i] = link_rule_spatial_4[i];
    }

    m->insert[0] = insert_rule_single_box_vertex_1;
    m->insert[1] = insert_rule_single_box_vertex_2;
    m->insert[2] = insert_rule_single_box_non_vertex;

    m->nsite = nsite;
    m->nbond = nbond;
    m->mhnspin = mhnspin;
    create_cmf(m->cmf,m->bond2weight,nbond);

    return m;
}

