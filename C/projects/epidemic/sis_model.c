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
int insert_rule_recover_3(int* state) {
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

/* graph name : sus_frozen
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
