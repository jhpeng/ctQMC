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

int link_rule_singlet_proj_1[] = {0,0,2,2,2,1,2,1};

int insert_rule_singlet_proj_1(int* state) {
    if(state[0]*state[1]==-1) {
        return 1;
    }
    return 0;
}

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
