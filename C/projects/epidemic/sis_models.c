#include <stdio.h>
#include <stdlib.h>

#include "dtype.h"

/* graph name : same_state
**     2    3
**     |    |
**     |____|
**     |    |
**     0    1
**
** density : 2*alpha
*/

int link_rule_same_state[] = {0,0,0,0,4,1,1,1};
int insert_rule_same_state(int* state) {
    if(state[0]==state[1]) {
        return 1;
    }
    return 0;
}

/* graph name : recover_1
**     2    3
**     o    o
**           
**     o    |
**     0    1
**
** density : gamma_j
*/

int link_rule_recover_1[] = {0,1,0,0,-3,1,-1,-1};
int insert_rule_recover_1(int* state) {
    if(state[0]==-1 && state[1]==-1) {
        return 1;
    }
    return 0;
}

/* graph name : recover_2
**     2    3
**     o    o
**           
**     |    o
**     0    1
**
** density : gamma_i
*/

int link_rule_recover_2[] = {0,1,1,1,1,-3,-1,-1};
int insert_rule_recover_2(int* state) {
    if(state[0]==-1 && state[1]==-1) {
        return 1;
    }
    return 0;
}

/* graph name : recover_3
**     2    3
**     x    o
**           
**     x    |
**     0    1
**
** density : gamma_j
*/

int link_rule_recover_3[] = {0,1,0,0,-3,1,-1,-1};
int insert_rule_recover_3(int* state) {
    if(state[0]==1 && state[1]==-1) {
        return 1;
    }
    return 0;
}

/* graph name : recover_4
**     2    3
**     o    x
**           
**     |    x
**     0    1
**
** density : gamma_i
*/

int link_rule_recover_4[] = {0,1,1,1,1,-3,-1,-1};
int insert_rule_recover_4(int* state) {
    if(state[0]==-1 && state[1]==1) {
        return 1;
    }
    return 0;
}

/* graph name : infect_1
**     2    3
**     x    |
**           
**     x    o
**     0    1
**
** density : alpha
*/

int link_rule_infect_1[] = {0,0,0,3,-3,-1,-1,1};
int insert_rule_infect_1(int* state) {
    if(state[0]==1 && state[1]==-1) {
        return 1;
    }
    return 0;
}

/* graph name : infect_2
**     2    3
**     |    x
**           
**     o    x
**     0    1
**
** density : alpha
*/

int link_rule_infect_2[] = {0,0,2,0,-3,-1,1,-1};
int insert_rule_infect_2(int* state) {
    if(state[0]==-1 && state[1]==1) {
        return 1;
    }
    return 0;
}

/* graph name : single_sie_recover_1
**      1
**      o
**            
**      |
**      0 
**
** density : 0.5*gamma
*/

int link_rule_single_site_recover_1[] = {0,1,1,-1};
int insert_rule_single_site_recover_1(int* state) {
    if(state[0]==-1) {
        return 1;
    }
    return 0;
}

/* graph name : single_sie_recover_2
**      1
**      |
**            
**      x
**      0 
**
** density : 0.5*gamma
*/

int link_rule_single_site_recover_2[] = {0,1,-1,1};
int insert_rule_single_site_recover_2(int* state) {
    if(state[0]==1) {
        return 1;
    }
    return 0;
}

/* graph name : susceptible_frozen
**      1
**      o
**        
**      o
**      0 
**
** density : 0.5*gamma
*/

int link_rule_susceptible_frozen[] = {0,0,-2,-1};
int insert_rule_susceptible_frozen(int* state) {
    if(state[0]==-1) {
        return 1;
    }
    return 0;
}

/* graph name : all_type_frozen
**      1
**      |
**      | 
**      |
**      0 
**
** apply on the boundary
*/

int link_rule_all_type_frozen[] = {0,0,-2,-1};
int insert_rule_all_type_frozen(int* state) {
    return 1;
}




static void create_cmf(double* cmf, double* weight, int length) {
    int i=0;

    cmf[0] = weight[i];
    for(i=1;i<length;i++) {
        cmf[i] = cmf[i-1]+weight[i];
    }
}

model* sis_model_uniform_infection(double alpha, int nnode, int nedge, int* edges) {
    int nsite = nnode;
    int nbond = 3*nedge+3*nnode;
    int mhnspin = 2;
    model* m = malloc_model(nsite,nbond+nnode,mhnspin);

    int n=0;
    for(int i_edge=0;i_edge<nedge;i_edge++) {
        int i = edges[2*i_edge+0];
        int j = edges[2*i_edge+1];

        m->bond2type[n]   = 0;
        m->bond2hNspin[n] = 2;
        m->bond2weight[n] = 2*alpha;
        m->sweight += 2*alpha;
        m->bond2index[n*mhnspin+0] = i;
        m->bond2index[n*mhnspin+1] = j;
        n++;
    }
    for(int i_edge=0;i_edge<nedge;i_edge++) {
        int i = edges[2*i_edge+0];
        int j = edges[2*i_edge+1];

        m->bond2type[n]   = 1;
        m->bond2hNspin[n] = 2;
        m->bond2weight[n] = alpha;
        m->sweight += alpha;
        m->bond2index[n*mhnspin+0] = i;
        m->bond2index[n*mhnspin+1] = j;
        n++;
    }
    for(int i_edge=0;i_edge<nedge;i_edge++) {
        int i = edges[2*i_edge+0];
        int j = edges[2*i_edge+1];

        m->bond2type[n]   = 2;
        m->bond2hNspin[n] = 2;
        m->bond2weight[n] = alpha;
        m->sweight += alpha;
        m->bond2index[n*mhnspin+0] = i;
        m->bond2index[n*mhnspin+1] = j;
        n++;
    }
    for(int i=0;i<nnode;i++) {

        m->bond2type[n]   = 3;
        m->bond2hNspin[n] = 1;
        m->bond2weight[n] = 0.5;
        m->sweight += 0.5;
        m->bond2index[n*mhnspin+0] = i;
        m->bond2index[n*mhnspin+1] = -1;
        n++;
    }
    for(int i=0;i<nnode;i++) {

        m->bond2type[n]   = 4;
        m->bond2hNspin[n] = 1;
        m->bond2weight[n] = 0.5;
        m->sweight += 0.5;
        m->bond2index[n*mhnspin+0] = i;
        m->bond2index[n*mhnspin+1] = -1;
        n++;
    }
    for(int i=0;i<nnode;i++) {

        m->bond2type[n]   = 5;
        m->bond2hNspin[n] = 1;
        m->bond2weight[n] = 1.0;
        m->sweight += 1.0;
        m->bond2index[n*mhnspin+0] = i;
        m->bond2index[n*mhnspin+1] = -1;
        n++;
    }

    /* type : 6
    ** apply only one the boundary
    */
    for(int i=0;i<nnode;i++) {

        m->bond2type[n]   = 6;
        m->bond2hNspin[n] = 1;
        m->bond2weight[n] = 0.0;
        m->sweight += 0.0;
        m->bond2index[n*mhnspin+0] = i;
        m->bond2index[n*mhnspin+1] = -1;
        n++;
    }

    for(int i=0;i<8;i++) {
        m->link[0*4*mhnspin+i] = link_rule_same_state[i];
        m->link[1*4*mhnspin+i] = link_rule_infect_1[i];
        m->link[2*4*mhnspin+i] = link_rule_infect_2[i];
    }
    for(int i=0;i<4;i++) {
        m->link[3*4*mhnspin+i] = link_rule_single_site_recover_1[i];
        m->link[4*4*mhnspin+i] = link_rule_single_site_recover_2[i];
        m->link[5*4*mhnspin+i] = link_rule_susceptible_frozen[i];
        m->link[6*4*mhnspin+i] = link_rule_all_type_frozen[i];
    }

    m->insert[0] = insert_rule_same_state;
    m->insert[1] = insert_rule_infect_1;
    m->insert[2] = insert_rule_infect_2;
    m->insert[3] = insert_rule_single_site_recover_1;
    m->insert[4] = insert_rule_single_site_recover_2;
    m->insert[5] = insert_rule_susceptible_frozen;
    m->insert[6] = insert_rule_all_type_frozen;

    m->nsite = nsite;
    m->nbond = nbond;
    m->mhnspin = mhnspin;
    create_cmf(m->cmf,m->bond2weight,nbond);

    return m;
}
