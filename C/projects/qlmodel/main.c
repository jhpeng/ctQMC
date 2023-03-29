#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <omp.h>

#include "dtype.h"
#include "models.h"
#include "update_omp.h"
#include "stats.h"
#include "union_find.h"

/**
 * @struct block
 * @brief Structure representing a block of data samples.
 *
 * @var block::samples
 * Pointer to an array of double values representing the samples.
 * @var block::nsample
 * Number of samples in the block.
 * @var block::index
 * Current index within the block.
 * @var block::cap
 * Maximum capacity of the block.
 */
typedef struct block{
    double* samples;
    int nsample;
    int index;
    int cap;
} block;

// parameter
int Lx;
int Ly;
double Lambda;
double Beta;
int Distance;
unsigned long int Seed;

// samples data
block* Nog0_block;
block* Nog1_block;
block* Nog2_block;
block* Order1_block;
block* Order2_block;

/**
 * @brief Allocates memory for a new block.
 *
 * @param cap The maximum capacity of the block.
 * @return A pointer to the newly allocated block.
 */
block* block_alloc(int cap) {
    block* data = (block*)malloc(sizeof(block));
    data->samples = (double*)malloc(sizeof(double)*cap);
    data->nsample = 0;
    data->index   = 0;
    data->cap     = cap;

    return data;
}

/**
 * @brief Frees memory of a block.
 *
 * @param data Pointer to the block to be freed.
 */
void block_free(block* data) {
    free(data->samples);
    free(data);
}

/**
 * @brief Appends a sample to the block.
 *
 * @param data Pointer to the block.
 * @param sample The sample to be appended.
 */
void append(block* data, double sample) {
    if(data->index == data->cap) 
        data->index = 0;
    
    data->samples[data->index] = sample;
    data->index++;

    if(data->nsample < data->cap) data->nsample++;
}

/**
 * @brief Computes the mean value of the samples in a block.
 *
 * @param data Pointer to the block.
 * @return The mean value of the samples in the block.
 */
double mean(block* data) {
    double sum=0;

    for(int i=0;i<(data->nsample);i++) {
        sum += data->samples[i];
    }

    return sum/(data->nsample);
}

/**
 * @brief Inserts a vertex into the world-line configuration to enforce Gauss law.
 *
 * @param w Pointer to the world_line_omp structure.
 * @param m Pointer to the model structure.
 * @param n Index of the vertex in the sequence.
 * @param bond Index of the bond.
 */
static void insert_gauss_law_graph(world_line_omp* w, model*m, int n, int bond) {
    vertex* sequence = w->sequenceB[w->nthread-1];
    if(w->flag[w->nthread-1])
        sequence = w->sequenceA[w->nthread-1];

    int u[4];
    vertex* v;
    v = &(sequence[n]);
    v->tau    = 2.0;
    v->bond   = bond;
    v->hNspin = 4;
    
    u[0] = m->bond2index[bond*(m->mhnspin)+0];
    u[1] = m->bond2index[bond*(m->mhnspin)+1];
    u[2] = m->bond2index[bond*(m->mhnspin)+2];
    u[3] = m->bond2index[bond*(m->mhnspin)+3];

    v->state[0] = w->istate[u[0]];
    v->state[1] = w->istate[u[1]];
    v->state[2] = w->istate[u[2]];
    v->state[3] = w->istate[u[3]];
    v->state[4] = w->istate[u[0]];
    v->state[5] = w->istate[u[1]];
    v->state[6] = w->istate[u[2]];
    v->state[7] = w->istate[u[3]];
}

/**
 * @brief Enforces Gauss law on a world-line configuration using parallel threads.
 *
 * @param w Pointer to the world_line_omp structure.
 * @param m Pointer to the model structure.
 * @param rng Pointer to an array of gsl_rng random number generator objects.
 */
static void gauss_law(world_line_omp* w, model* m, gsl_rng** rng) {
    int lsize = Lx*Ly;
    int len = (w->len[w->nthread-1]) + lsize;
    realloc_world_line_omp_vertex(w, len, w->nthread-1);

    omp_set_num_threads(w->nthread);
    #pragma omp parallel
    {
        int i_thread = omp_get_thread_num();
        int start = lsize*i_thread/(w->nthread);
        int end = lsize*(i_thread+1)/(w->nthread);

        int u1,u2,u3,u4,bond,n;
        int state[4];
        int charge;

        for(int i_site=start; i_site<end; i_site++) {
            n = w->len[w->nthread-1]+i_site;
            bond = 3*lsize+4*i_site+3;
            u1 = m->bond2index[bond*(m->mhnspin)+0];
            u2 = m->bond2index[bond*(m->mhnspin)+1];
            u3 = m->bond2index[bond*(m->mhnspin)+2];
            u4 = m->bond2index[bond*(m->mhnspin)+3];

            state[0] = w->istate[u1];
            state[1] = w->istate[u2];
            state[2] = w->istate[u3]*(-1);
            state[3] = w->istate[u4]*(-1);

            charge = state[0]+state[1]+state[2]+state[3];

            if(charge==0) {
                if(state[0]==state[1]) {
                    if(gsl_rng_uniform_pos(rng[i_thread])<0.5){
                        bond = 3*lsize+4*i_site+0;
                    } else {
                        bond = 3*lsize+4*i_site+2;
                    }
                } else if(state[0]==state[2]) {
                    if(gsl_rng_uniform_pos(rng[i_thread])<0.5){
                        bond = 3*lsize+4*i_site+0;
                    } else {
                        bond = 3*lsize+4*i_site+1;
                    }
                }else if(state[0]==state[3]) {
                    if(gsl_rng_uniform_pos(rng[i_thread])<0.5){
                        bond = 3*lsize+4*i_site+1;
                    } else {
                        bond = 3*lsize+4*i_site+2;
                    }
                }
            } else {
                bond = 3*lsize+4*i_site+3;
            }

            insert_gauss_law_graph(w, m, n, bond);
        }
    }

    w->len[w->nthread-1] += lsize;
}

/**
 * @brief Prints the charge configuration of the system.
 *
 * This function calculates and prints the charge configuration of the system, based on the
 * world line data. The charges are calculated for each site of the lattice and displayed
 * according to the following format:
 *  0: no charge
 *  -: negative charge (charge = -2)
 *  +: positive charge (charge = 2)
 *  2: double charge (charge = -4 or 4)
 *
 * @param w Pointer to the world_line_omp structure containing the world line data.
 */
void print_charge(world_line_omp* w) {
    int state[4],charge;
    for(int y=0;y<Ly;y++) {
        for(int x=0;x<Lx;x++) {
            int u1 = x+y*Lx;
            int u2 = x+y*Lx+Lx*Ly;
            int u3 = (x+Lx-1)%Lx+y*Lx;
            int u4 = x+((y+Ly-1)%Ly)*Lx+Lx*Ly;

            state[0] = w->istate[u1];
            state[1] = w->istate[u2];
            state[2] = w->istate[u3]*(-1);
            state[3] = w->istate[u4]*(-1);

            charge = state[0]+state[1]+state[2]+state[3];

            if(charge==0) printf("%d ",charge);
            else if(charge==-2) printf("- ");
            else if(charge==2) printf("+ ");
            else if(charge==-4 || charge==4) printf("2 ");
        }
        printf("\n");
    }
}

/**
 * @brief Initializes the world line with no charge configuration.
 *
 * This function initializes the world line with a configuration that has no charge.
 * The function takes two arguments, wx and wy, which are used to adjust the initial state
 * of the world line.
 *
 * @param w Pointer to the world_line_omp structure to be initialized.
 * @param wx Width of the region in the x direction to be set to 1.
 * @param wy Width of the region in the y direction to be set to 1.
 */
void initial_state_no_charge(world_line_omp* w, int wx, int wy){
    for(int y=0;y<Ly;y++) {
        for(int x=0;x<Lx;x++) {
            w->istate[x+y*Lx]         = 1-(x+y)%2*2;
            w->istate[x+y*Lx + Lx*Ly] = (x+y)%2*2-1;
        }
    }

    for(int y=0;y<wx;y++) {
        for(int x=0;x<Lx;x++) {
            w->istate[x+y*Lx] = 1;
        }
    }

    for(int y=0;y<Ly;y++) {
        for(int x=0;x<wy;x++) {
            w->istate[x+y*Lx+Lx*Ly] = 1;
        }
    }

    for(int j=1;j<w->nthread;j++) {
        for(int i=0;i<(w->nsite);i++) {
            w->istate[i+(w->nsite)*j] = w->istate[i];
        }
    }
}

/**
 * @brief Initializes the world line with charge 1 configuration.
 *
 * This function initializes the world line with a configuration that has a charge of 1.
 * The function takes one argument, distance, which is used to adjust the initial state
 * of the world line.
 *
 * @param w Pointer to the world_line_omp structure to be initialized.
 * @param distance Distance between charges in the initial state.
 */
void initial_state_charge_1(world_line_omp* w, int distance) {
    for(int y=0;y<Ly;y++) {
        for(int x=0;x<Lx;x++) {
            w->istate[x+y*Lx]         = 1-(y%2)*2;
            w->istate[x+y*Lx + Lx*Ly] = (x%2)*2-1;
        }
    }

    for(int x=(Lx-distance)/2;x<(Lx+distance)/2;x++) {
        int y = Ly/2;
        w->istate[x+y*Lx] *= -1;
    }

    for(int j=1;j<w->nthread;j++) {
        for(int i=0;i<(w->nsite);i++) {
            w->istate[i+(w->nsite)*j] = w->istate[i];
        }
    }
}

/**
 * @brief Initializes the world line with charge 1 configuration in diagonal direction.
 *
 * This function initializes the world line with a configuration that has a charge of 1
 * in the diagonal direction. The function takes one argument, distance, which is used
 * to adjust the initial state of the world line.
 *
 * @param w Pointer to the world_line_omp structure to be initialized.
 * @param distance Distance between charges in the initial state.
 */
void initial_state_charge_1_diag(world_line_omp* w, int distance) {
    for(int y=0;y<Ly;y++) {
        for(int x=0;x<Lx;x++) {
            w->istate[x+y*Lx]         = 1-(y%2)*2;
            w->istate[x+y*Lx + Lx*Ly] = (x%2)*2-1;
        }
    }

    for(int x=(Lx-distance)/2;x<(Lx+distance)/2;x++) {
        int y = (Ly-distance)/2;
        w->istate[x+y*Lx] *= -1;
    }

    for(int y=(Lx-distance)/2;y<(Lx+distance)/2;y++) {
        int x = (Ly-distance)/2;
        w->istate[x+y*Lx+Lx*Ly] *= -1;
    }

    for(int j=1;j<w->nthread;j++) {
        for(int i=0;i<(w->nsite);i++) {
            w->istate[i+(w->nsite)*j] = w->istate[i];
        }
    }
}

/**
 * @brief Checks if the input state is a reference configuration.
 *
 * This function checks if the input state is a reference configuration.
 *
 * @param state Pointer to an integer array representing the input state.
 * @return 1 if the input state is a reference configuration, otherwise 0.
 */
static int reference_conf(int* state){
    if(state[0]*state[1]==1){
        if(state[2]*state[3]==1){
            if(state[1]*state[2]==-1){
                return 1;
            }
        }
    }
    return 0;
}

/**
 * @brief Counts the local graph number and updates the corresponding arrays.
 *
 * This function counts the local graph number and updates the corresponding arrays
 * nog0_local_number, nog1_local_number, and nog2_local_number.
 *
 * @param w Pointer to the world_line_omp structure containing the world line data.
 * @param m Pointer to the model structure containing the model data.
 */
int count_local_number=0;
unsigned long int* nog0_local_number;
unsigned long int* nog1_local_number;
unsigned long int* nog2_local_number;

void count_local_graph_number(world_line_omp* w, model* m) {
    int nthread = w->nthread;
    int lsize = Lx*Ly;

    if(count_local_number==0) {
        nog0_local_number = (unsigned long int*)malloc(sizeof(unsigned long int)*lsize*nthread);
        nog1_local_number = (unsigned long int*)malloc(sizeof(unsigned long int)*lsize*nthread);
        nog2_local_number = (unsigned long int*)malloc(sizeof(unsigned long int)*lsize*nthread);

        for(int i=0;i<(nthread*lsize);i++) {
            nog0_local_number[i] = 0;
            nog1_local_number[i] = 0;
            nog2_local_number[i] = 0;
        }
    }

    omp_set_num_threads(nthread);
    #pragma omp parallel
    {
        int i_thread = omp_get_thread_num();
        int type,bond;

        vertex* sequence = w->sequenceB[i_thread];
        if(w->flag[i_thread])
            sequence = w->sequenceA[i_thread];

        for(int i=0;i<(w->len[i_thread]);i++) {
            bond = (sequence[i]).bond;
            type = m->bond2type[bond];

            if(type==0) {
                nog0_local_number[i_thread*lsize+(bond%lsize)]++;
            } else if(type==1) {
                nog1_local_number[i_thread*lsize+(bond%lsize)]++;
            } else if(type==2) {
                nog2_local_number[i_thread*lsize+(bond%lsize)]++;
            }
        }
    }

    count_local_number++;
}

/**
 * @brief Saves the local energy of the system to text files.
 *
 * This function saves the local energy of the system to three text files named
 * "nog0", "nog1", and "nog2", corresponding to the three different types of local
 * graph numbers. The function uses the global parameters Lx, Ly, Lambda, Beta,
 * Distance, and Seed to create the file names.
 *
 * @param w Pointer to the world_line_omp structure containing the world line data.
 */
void save_local_energy(world_line_omp* w) {
    int nthread = w->nthread;
    int lsize = Lx*Ly;

    for(int i=0;i<lsize;i++) {
        for(int i_thread=1;i_thread<nthread;i_thread++) {
            nog0_local_number[i] += nog0_local_number[i_thread*lsize+i];
            nog1_local_number[i] += nog1_local_number[i_thread*lsize+i];
            nog2_local_number[i] += nog2_local_number[i_thread*lsize+i];
        }
    }

    char nog0_file_name[128];
    char nog1_file_name[128];
    char nog2_file_name[128];

    sprintf(nog0_file_name,"map/nog0_Lx_%d_Ly_%d_lambda_%.4f_beta_%.2f_distance_%d_seed_%ld.txt",Lx,Ly,Lambda,Beta,Distance,Seed);
    sprintf(nog1_file_name,"map/nog1_Lx_%d_Ly_%d_lambda_%.4f_beta_%.2f_distance_%d_seed_%ld.txt",Lx,Ly,Lambda,Beta,Distance,Seed);
    sprintf(nog2_file_name,"map/nog2_Lx_%d_Ly_%d_lambda_%.4f_beta_%.2f_distance_%d_seed_%ld.txt",Lx,Ly,Lambda,Beta,Distance,Seed);

    FILE* nog0_file = fopen(nog0_file_name,"w");
    FILE* nog1_file = fopen(nog1_file_name,"w");
    FILE* nog2_file = fopen(nog2_file_name,"w");

    for(int y=0;y<Ly;y++) {
        for(int x=0;x<Lx;x++) {
            double local_energy;
            local_energy = ((double)nog0_local_number[x+y*Lx])/Beta/count_local_number;
            fprintf(nog0_file,"%.16e ",local_energy);

            local_energy = ((double)nog1_local_number[x+y*Lx])/Beta/count_local_number;
            fprintf(nog1_file,"%.16e ",local_energy);

            local_energy = ((double)nog2_local_number[x+y*Lx])/Beta/count_local_number;
            fprintf(nog2_file,"%.16e ",local_energy);
        }
        fprintf(nog0_file,"\n");
        fprintf(nog1_file,"\n");
        fprintf(nog2_file,"\n");
    }

    free(nog0_local_number);
    free(nog1_local_number);
    free(nog2_local_number);
    count_local_number = 0;

    fclose(nog0_file);
    fclose(nog1_file);
    fclose(nog2_file);
}

/**
 * @brief Counts the number of reference configurations in the given world line.
 *
 * @param[in] w Pointer to the world_line_omp structure.
 * @return The number of reference configurations.
 */
int count_reference_conf(world_line_omp* w) {
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

            if(reference_conf(state)) count++;
        }
    }

    return count;
}

/**
 * @brief Counts the graph numbers for each type in the given world line and model.
 *
 * @param[in] w Pointer to the world_line_omp structure.
 * @param[in] m Pointer to the model structure.
 * @param[out] nog Array of size 3 to store the graph numbers for each type.
 */
void count_graph_number(world_line_omp* w, model* m, int* nog) {
    int nthread = w->nthread;
    int nog0[nthread];
    int nog1[nthread];
    int nog2[nthread];

    omp_set_num_threads(nthread);
    #pragma omp parallel
    {
        int i_thread = omp_get_thread_num();
        int type,bond;

        nog0[i_thread] = 0;
        nog1[i_thread] = 0;
        nog2[i_thread] = 0;

        vertex* sequence = w->sequenceB[i_thread];
        if(w->flag[i_thread])
            sequence = w->sequenceA[i_thread];

        for(int i=0;i<(w->len[i_thread]);i++) {
            bond = (sequence[i]).bond;
            type = m->bond2type[bond];

            if(type==0) nog0[i_thread]++;
            else if(type==1) nog1[i_thread]++;
            else if(type==2) nog2[i_thread]++;
        }
    }

    nog[0]=0;
    nog[1]=0;
    nog[2]=0;
    for(int i_thread=0;i_thread<nthread;i_thread++) {
        nog[0] += nog0[i_thread];
        nog[1] += nog1[i_thread];
        nog[2] += nog2[i_thread];
    }
}

/**
 * @brief Measures the order of the given world line and stores the results in the obs array.
 *
 * @param[in] w Pointer to the world_line_omp structure.
 * @param[out] obs Array of size 2 to store the order measurement results.
 */
int* order1_structure_factor = NULL;
int* order2_structure_factor = NULL;
void measure_order(world_line_omp* w, double* obs) {
    if(order1_structure_factor==NULL) {
        order1_structure_factor = (int*)malloc(sizeof(int)*Lx*Ly*2);
        order2_structure_factor = (int*)malloc(sizeof(int)*Lx*Ly*2);

        for(int y=0;y<Ly;y++) {
            for(int x=0;x<Lx;x++) {
                order1_structure_factor[x+y*Lx]         = 1-(y%2)*2;
                order1_structure_factor[x+y*Lx + Lx*Ly] = 0;
            }
        }

        for(int y=0;y<Ly;y++) {
            for(int x=0;x<Lx;x++) {
                order2_structure_factor[x+y*Lx]         = 1-(x+y)%2*2;
                order2_structure_factor[x+y*Lx + Lx*Ly] = (x+y)%2*2-1;
            }
        }
    }

    int nthread = w->nthread;
    int order1_temp[nthread];
    int order2_temp[nthread];

    omp_set_num_threads(nthread);
    #pragma omp parallel
    {
        int i_thread = omp_get_thread_num();
        
        int start = 2*Lx*Ly*i_thread/nthread;
        int end   = 2*Lx*Ly*(i_thread+1)/nthread;

        order1_temp[i_thread] = 0;
        order2_temp[i_thread] = 0;

        for(int i=start;i<end;i++) {
            order1_temp[i_thread] += w->istate[i]*order1_structure_factor[i];
            order2_temp[i_thread] += w->istate[i]*order2_structure_factor[i];
        }
    }

    obs[0] = 0;
    obs[1] = 0;
    for(int i_thread=0;i_thread<nthread;i_thread++) {
        obs[0] += (double)order1_temp[i_thread];
        obs[1] += (double)order2_temp[i_thread];
    }
    obs[0] = fabs(obs[0])/(Lx*Ly);
    obs[1] = fabs(obs[1])/(2*Lx*Ly);
}

/**
 * @brief Performs measurement of the given world line and model.
 *
 * @param[in] w Pointer to the world_line_omp structure.
 * @param[in] m Pointer to the model structure.
 */
void measurement(world_line_omp* w, model* m) {
    if(0) {
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

        int n_ref_conf = count_reference_conf(w);
        printf("nref : %d \n",n_ref_conf);

        double energy=0;
        for(int j=0;j<w->nthread;j++)
            energy+=w->len[j];
        energy = -(energy-Lx*Ly)/Beta;

        printf("energy : %e \n",energy);
    }

    int nog[3];
    count_graph_number(w,m,nog);
    append(Nog0_block,((double)nog[0])/Beta);
    append(Nog1_block,((double)nog[1])/Beta);
    append(Nog2_block,((double)nog[2])/Beta);

    double obs_order[2];
    measure_order(w,obs_order);
    append(Order1_block,obs_order[0]);
    append(Order2_block,obs_order[1]);


}

/**
 * @brief Saves the measurement data to a file.
 */
void save_data() {
    char data_file_name[128];
    sprintf(data_file_name,"data/data_Lx_%d_Ly_%d_lambda_%.4f_beta_%.2f_distance_%d_seed_%ld.txt",Lx,Ly,Lambda,Beta,Distance,Seed);
    FILE* data_file = fopen(data_file_name,"a");


    fprintf(data_file,"%.16e ",mean(Nog0_block));
    fprintf(data_file,"%.16e ",mean(Nog1_block));
    fprintf(data_file,"%.16e ",mean(Nog2_block));
    fprintf(data_file,"%.16e ",mean(Order1_block));
    fprintf(data_file,"%.16e ",mean(Order2_block));
    fprintf(data_file,"\n");

    fclose(data_file);
}

int main(int argc, char** argv) {
    Lx = 128;
    Ly = 128;
    Lambda = 1.0;
    Beta = 128.0;
    Distance = 0;
    Seed = 47583;

    int nthread = 6;
    int thermal = 100000;
    int nsample = 1000;
    int nblock  = 1000;

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
    //initial_state_no_charge(w,0,0);
    //initial_state_charge_1(w,Distance);
    initial_state_charge_1_diag(w,Distance);

    // Settinh the samples data
    Nog0_block = block_alloc(nsample);
    Nog1_block = block_alloc(nsample);
    Nog2_block = block_alloc(nsample);
    Order1_block = block_alloc(nsample);
    Order2_block = block_alloc(nsample);

    // Thermalization
    int ntimes = 100;
    double* times = (double*)malloc(sizeof(double)*ntimes);
    for(int i=0;i<thermal;i++) {
        double start = omp_get_wtime();
        remove_vertices_omp(w);
        double t1 = omp_get_wtime();
        insert_vertices_omp(w,m,rng);
        gauss_law(w,m,rng);
        double t2 = omp_get_wtime();
        clustering_inner_omp(w,m);
        double t3 = omp_get_wtime();
        clustering_crossing_omp(w);
        double t4 = omp_get_wtime();
        flip_cluster_omp(w,rng);
        double end = omp_get_wtime();
        //check_world_line_omp_configuration(w,m);

        measurement(w,m);

        times[i%ntimes] = end-start;

        //if((i+1)%ntimes==0) {
        if(1) {

        int noo=0;
        for(int j=0;j<nthread;j++)
            noo+=w->len[j];


        printf("-----------------------------------------\n");
        printf("Lx = %d | Ly = %d | lambda=%.4f | beta=%.1f \n",Lx,Ly,Lambda,Beta);
        printf("thremal:%d | Noo:%d | nthread=%d ",i+1,noo,nthread);
        if(i>50) {
            double mtime=0;
            for(int j=0;j<ntimes;j++) mtime+=times[j];
            mtime = mtime/ntimes;
            printf("| time:%f \n",mtime);
        } else {
            printf("\n");
        }
        printf("remove_vertices_omp     : %f  \n",(t1-start)/(end-start));
        printf("insert_vertices_omp     : %f  \n",(t2-t1)/(end-start));
        printf("clustering_inner_omp    : %f  \n",(t3-t2)/(end-start));
        printf("clustering_crossing_omp : %f  \n",(t4-t3)/(end-start));
        printf("flip_cluster_omp        : %f  \n",(end-t4)/(end-start));
        printf("order1 : %e \n",mean(Order1_block));
        printf("order2 : %e \n",mean(Order2_block));

        //print_charge(w);
        }

    }

    for(int i_block=0;i_block<nblock;i_block++) {
        for(int i=0;i<nsample;i++) {
            remove_vertices_omp(w);
            insert_vertices_omp(w,m,rng);
            gauss_law(w,m,rng);
            clustering_inner_omp(w,m);
            clustering_crossing_omp(w);
            flip_cluster_omp(w,rng);

            measurement(w,m);
            count_local_graph_number(w,m);
            //print_charge(w);
        }

        int noo=0;
        for(int j=0;j<nthread;j++)
            noo+=w->len[j];

        printf("-----------------------------------------\n");
        printf("Lx = %d | Ly = %d | lambda=%.4f | beta=%.1f \n",Lx,Ly,Lambda,Beta);
        printf("nblock:%d | Noo:%d | nthread=%d \n",i_block+1,noo,nthread);

        printf("nog0/beta : %e \n",mean(Nog0_block));
        printf("nog1/beta : %e \n",mean(Nog1_block));
        printf("nog2/beta : %e \n",mean(Nog2_block));
        printf("order1 : %e \n",mean(Order1_block));
        printf("order2 : %e \n",mean(Order2_block));

        //save_data();
    }

    //save_local_energy(w);

    //free data
    for(int i_thread=0;i_thread<nthread;i_thread++)
        gsl_rng_free(rng[i_thread]);
    free(times);
    free_model(m);
    free_world_line_omp(w);

    block_free(Nog0_block);
    block_free(Nog1_block);
    block_free(Nog2_block);
    block_free(Order1_block);
    block_free(Order2_block);

}
