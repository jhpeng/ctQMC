#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "sis_models.h"
#include "update.h"
#include "networks.h"
#include "estimator.h"

// Returns the initial number of infected nodes in the world line
int ninfected_initial_state(world_line* w) {
    int nnode=w->nsite;  // Number of nodes in the world line
    int ninfected=0;     // Number of infected nodes

    // Iterate over all nodes in the world line
    for(int i=0;i<nnode;i++) {
        ninfected+=w->istate[i];  // Add 1 to ninfected for each infected node
    }

    // Return the average of ninfected and nnode, rounded up
    return (ninfected+nnode)/2;
}

// Returns the final number of infected nodes in the world line
int ninfected_final_state(world_line* w) {
    int nnode=w->nsite;  // Number of nodes in the world line
    int ninfected=0;     // Number of infected nodes

    // Iterate over all nodes in the world line
    for(int i=0;i<nnode;i++) {
        ninfected+=w->pstate[i];  // Add 1 to ninfected for each infected node
    }

    // Return the average of ninfected and nnode, rounded up
    return (ninfected+nnode)/2;
}


// Define a global variable to store the frozen list for boundary condition type 1
int* boundary_condition_frozen_list=NULL;

void boundary_condition_initial_state(world_line* w, model* m, int type, gsl_rng* rng) {
    // Get the number of nodes and bonds in the model
    int nnode = m->nsite;
    int nbond = m->nbond;

    // Calculate the length of the sequence after adding nodes
    int length = nnode+(w->nvertices);

    // Resize the world line sequence to fit the new nodes
    realloc_world_line(w,length);

    // Pointers to the two sequences of vertices in the world line
    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;

    // Swap the pointers if the flag in the world line is true
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    // Allocate memory for the frozen list if it hasn't been allocated yet
    if(boundary_condition_frozen_list==NULL) {
        boundary_condition_frozen_list = (int*)malloc(sizeof(int)*nnode);
    }

    // Counter for the number of vertices in the sequence
    int n=0;

    // Assign initial states to the vertices based on the boundary condition type
    if(type==0) {
        for(int i=0;i<nnode;i++) {
            (sequence2[n]).tau      = 0.0;
            (sequence2[n]).bond     = nbond+i;
            (sequence2[n]).hNspin   = 1;
            (sequence2[n]).state[0] = w->istate[i];
            (sequence2[n]).state[1] = w->istate[i];
            n++;
        }
    } else if(type==1) {
        int ninfected=0;
        int i_node=-1;
        for(int i=0;i<nnode;i++) {
            if(w->istate[i]==-1) {
                boundary_condition_frozen_list[i]=1;
            } else {
                boundary_condition_frozen_list[i]=0;
                i_node=i;
                ninfected++;
            }
        }
        if(ninfected==1) {
            int check=1;
            while(check) {
                int j = nearest_nb_random_assign(i_node,rng);
                if(boundary_condition_frozen_list[j]) {
                    boundary_condition_frozen_list[j]=0;
                    check=0;
                }
            }
        } else if(ninfected==0) {
            for(int i=0;i<nnode;i++) 
                boundary_condition_frozen_list[i]=0;
        }
        for(int i=0;i<nnode;i++) {
            if(boundary_condition_frozen_list[i]) {
                (sequence2[n]).tau      = 0.0;
                (sequence2[n]).bond     = nbond+i;
                (sequence2[n]).hNspin   = 1;
                (sequence2[n]).state[0] = w->istate[i];
                (sequence2[n]).state[1] = w->istate[i];
                n++;
            }
        }
    }


    for(int i=0;i<(w->nvertices);i++) {
        copy_vertex(&(sequence2[n]),&(sequence1[i]));
        n++;
    }

    w->nvertices = n;
    w->flag = !(w->flag);
}

void boundary_condition_final_state(world_line* w, model* m, double p, int type, gsl_rng* rng) {
    // get the number of nodes and bonds, and calculate the required length
    int nnode = m->nsite;
    int nbond = m->nbond;
    int length = nnode+(w->nvertices);
    // allocate or reallocate memory for the world line
    realloc_world_line(w,length);

    // get the appropriate sequence depending on the flag
    vertex* sequence = w->sequenceB;
    if(w->flag) {
        sequence = w->sequenceA;
    }

    // get the pointer to the pstate array, and the number of vertices
    int* pstate = w->pstate;
    int n=w->nvertices;

    // set the appropriate final state based on the type of boundary condition
    if(type==0) {
        // set all nodes to their initial state with tau=1
        for(int i=0;i<nnode;i++) {
            (sequence[n]).tau      = 1.0;
            (sequence[n]).bond     = nbond+i;
            (sequence[n]).hNspin   = 1;
            (sequence[n]).state[0] = pstate[i];
            (sequence[n]).state[1] = pstate[i];
            n++;
        }
    } else if(type==1) {
        // set boundary nodes based on the infection state of the system and the probability p
        int inf=0;
        for(int i=0;i<nnode;i++) inf += (pstate[i]+1)/2;

        double pdis = 1.0;
        if(inf!=0) pdis=(p*nnode)/inf;
        if(inf>p*nnode) pdis=0;
        
        for(int i=0;i<nnode;i++) {
            if(pstate[i]==1 && (gsl_rng_uniform_pos(rng)<pdis)) {
                (sequence[n]).tau      = 1.0;
                (sequence[n]).bond     = nbond+i;
                (sequence[n]).hNspin   = 1;
                (sequence[n]).state[0] = pstate[i];
                (sequence[n]).state[1] = pstate[i];
                n++;
            }
        }
    } else if(type==2) {
        // set frozen nodes to their final state with tau=1
        for(int i=0;i<nnode;i++) {
            if(pstate[i]==-1) {
                (sequence[n]).tau      = 1.0;
                (sequence[n]).bond     = nbond+i;
                (sequence[n]).hNspin   = 1;
                (sequence[n]).state[0] = pstate[i];
                (sequence[n]).state[1] = pstate[i];
                n++;
            }
        }
    }

    // update the number of vertices in the world line
    w->nvertices = n;
}

static void print_state(int* state, int nnode) {
    // This function prints the state of each node in the world_line.
    // It takes as input an array of node states and the number of nodes.
    for(int i=0;i<nnode;i++) {
        printf("%d ",(state[i]+1)/2);
    }
    printf("\n");
}

static void save_state(FILE* file, int* state, int nnode) {
    // This function writes the state of each node in the world_line to a file.
    // It takes as input a file pointer, an array of node states, and the number of nodes.
    for(int i=0;i<nnode;i++) {
        fprintf(file,"%d ",(state[i]+1)/2);
    }
    fprintf(file,"\n");
}

void show_configuration(world_line* w, model* m, double* time_list, int ntime) {
    // This function displays the configuration of the world_line at specified times.
    // It takes as input a pointer to the world_line, a pointer to the model, an array of times,
    // and the number of times to display.
    int* pstate = w->pstate;
    int nnode = w->nsite;
    int mhnspin = m->mhnspin;

    // Copy the initial state of the world_line nodes to pstate.
    for(int i=0;i<nnode;i++) pstate[i] = w->istate[i];

    // Get the vertex sequence to traverse.
    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    int i=0; // Index of current time to display
    int index, i_node; 
    vertex* v;
    for(int n=0;n<(w->nvertices) && i<ntime;n++) {
        v = &(sequence[n]);
        if(time_list[i]<(v->tau)) { // Check if it's time to display the state
            print_state(pstate,nnode);
            i++;
        }

        // Update the state of nodes according to vertex information.
        for(i_node=0;i_node<(v->hNspin);i_node++) {
            index = m->bond2index[(v->bond)*mhnspin+i_node];
            pstate[index] = v->state[(v->hNspin)+i_node];
        }
    }

    // If there are remaining times to display, display the last state.
    for(;i<ntime;i++) {
        print_state(pstate,nnode);
    }
}


void save_configuration(FILE* file, world_line* w, model* m, double* time_list, int ntime) {
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
            save_state(file,pstate,nnode);
            i++;
        }

        for(i_node=0;i_node<(v->hNspin);i_node++) {
            index = m->bond2index[(v->bond)*mhnspin+i_node];
            pstate[index] = v->state[(v->hNspin)+i_node];
        }
    }
    for(;i<ntime;i++) {
        save_state(file,pstate,nnode);
    }
}

static time_t start_time,end_time;
static unsigned long int measurement_count=0;
static double* infected_ratio=NULL;
static double* infected_time=NULL;
static double total_infected_time_ave=0;
static double ninfection_ave=0;
static double nrecover_ave=0;
static double ntrial_ave=0;
void measurement(world_line* w, model* m, double* time_list, int ntime, int block_size) {
    if(infected_ratio==NULL) {
        infected_time = (double*)malloc(sizeof(double)*(w->nsite));
        infected_ratio = (double*)malloc(sizeof(double)*ntime);
        for(int i=0;i<(w->nsite);i++) infected_time[i]=0;
        for(int i=0;i<ntime;i++) infected_ratio[i]=0;

        sequence_malloc(block_size,3);
        start_time = clock();
    }
    int* pstate = w->pstate;
    int nnode = w->nsite;
    int mhnspin = m->mhnspin;

    for(int i=0;i<nnode;i++) {
        pstate[i] = w->istate[i];
        infected_time[i] = 0;
    }

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    double total_infected_time=0;

    int i=0;
    int index, i_node;
    vertex* v;
    double tau_p=0;
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
            if(pstate[index]==1) {
                total_infected_time += ((v->tau)-infected_time[index]);
            }
            pstate[index] = v->state[(v->hNspin)+i_node];
            if(pstate[index]==1) {
                infected_time[index] = v->tau;
            }
        }
        if(tau_p>(v->tau)) printf("tau_p > tau!\n");
        else if((v->tau)<0) printf("tau < 0!\n");
        else if((v->tau)>1) printf("tau > 1!\n");
        tau_p = v->tau;
    }
    for(i_node=0;i_node<nnode;i_node++) {
        if(pstate[i_node]==1)
            total_infected_time += (1.0-infected_time[i_node]);
    }

    for(;i<ntime;i++) {
        double ir=0;
        for(int i_node=0;i_node<nnode;i_node++) {
            ir+=0.5*(pstate[i_node]+1);
        }
        ir = ir/nnode;
        infected_ratio[i] += ir;
    }

    // collecting the obeservable
    total_infected_time_ave += total_infected_time;
    ninfection_ave += ninfection_value();
    nrecover_ave  += nrecover_value();
    
    double samples[3];
    samples[0] = ninfection_value();
    samples[1] = nrecover_value();
    samples[2] = total_infected_time*(w->beta);
    sequence_append(samples);

    measurement_count++;

    if(measurement_count==block_size) {
        FILE* file_conf = fopen("conf.txt","a");
        FILE* file_t = fopen("times.txt","w");
        FILE* file_s = fopen("series.txt","a");
        FILE* file_g = fopen("global.txt","a");
        printf("------------------------------\n");
        printf(" t    |    I/N\n");
        for(i=0;i<ntime;i++) {
            infected_ratio[i] = infected_ratio[i]/block_size;
            printf("%.4lf  %.12lf\n",time_list[i]*w->beta,infected_ratio[i]);
            fprintf(file_t,"%.4lf ",time_list[i]*w->beta);
            fprintf(file_s,"%.12e ",infected_ratio[i]);

            infected_ratio[i] = 0;
        }
        fprintf(file_t,"\n");
        fprintf(file_s,"\n");
        measurement_count=0;

        nrecover_ave = nrecover_ave/block_size;
        ninfection_ave = ninfection_ave/block_size;
        ntrial_ave = ntrial_ave/block_size;
        total_infected_time_ave = total_infected_time_ave/block_size*(w->beta);
        fprintf(file_g,"%.12e %.12e %.12e %.12e\n",ninfection_ave,nrecover_ave,total_infected_time_ave,ntrial_ave);

        printf("total infected time = %.12e\n",total_infected_time_ave);
        printf("average # of trial  = %.12e\n",ntrial_ave);

        save_configuration(file_conf,w,m,time_list,ntime);

        total_infected_time_ave=0;
        nrecover_ave=0;
        ninfection_ave=0;
        ntrial_ave=0;

        fclose(file_conf);
        fclose(file_t);
        fclose(file_s);
        fclose(file_g);

        end_time = clock();
        printf("time for this block = %.2lf(sec)\n",(double)(end_time-start_time)/CLOCKS_PER_SEC);
        start_time = clock();
    }
}

int main(int argc, char** argv) {
    char filename[128] = "/hpc/home/jp549/src/ctQMC/C/projects/epidemic/network/test.edgelist";
    double alpha=atof(argv[1]);
    double gamma=atof(argv[2]);
    double T = atof(argv[3]);
    int nif=atoi(argv[4]);
    int running_mode=atoi(argv[5]);
    int block_size=atoi(argv[6]);
    int nblock = atoi(argv[7]);
    int thermal = atoi(argv[8]);
    int nskip = atoi(argv[9]);
    int nsweep  = nblock*block_size;
    unsigned long int seed=atoi(argv[10]);

    int nnode;
    int nedge;
    int* edges = read_edgelist(filename,&nnode,&nedge);
    double pnif = ((double)nif)/nnode;


    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);


    model* m = sis_model_uniform_infection(alpha,gamma,nnode,nedge,edges);
    world_line* w = malloc_world_line(1024,2*(m->mhnspin),m->nsite);

    // setup running mode and initial state
    // running mode : 0
    //      initial state - fixing patient0
    //      final state   - # of infection > nif
    //
    // running mode : 1
    //      initial state - moving patient0
    //      final state   - # of infection > nif
    //
    // running mode : 2 
    //      initial state - all get infected
    //      final state   - all get recovery
    //

    int initial_condition_type=0;
    int final_condition_type=0;
    int nocheck_for_measurement=0;

    if(running_mode==0) {
        for(int i=0;i<(w->nsite);i++) w->istate[i] = -1;
        w->istate[nearest_nb_arg_max_degree()]=1;
        //w->istate[71]=1;

        initial_condition_type=0;
        final_condition_type=1;
        nocheck_for_measurement=0;
    } else if(running_mode==1) {
        for(int i=0;i<(w->nsite);i++) w->istate[i] = -1;
        w->istate[nearest_nb_arg_max_degree()]=1;
        //w->istate[71]=1;

        initial_condition_type=1;
        final_condition_type=1;
        nocheck_for_measurement=0;
    } else if(running_mode==2) {
        for(int i=0;i<(w->nsite);i++) w->istate[i] = 1;

        initial_condition_type=0;
        final_condition_type=2;
        nocheck_for_measurement=1;
    } else {
        printf("There is no such running mode %d!\n",running_mode);
        exit(1);
    }

    // thermalization
    w->beta = T;
    for(int i=0;i<thermal;i++) {
        remove_vertices(w);
        swapping_graphs(w,m,rng);
        insert_vertices(w,m,rng);
        boundary_condition_initial_state(w,m,initial_condition_type,rng);
        boundary_condition_final_state(w,m,pnif,final_condition_type,rng);
        clustering(w,m);

        //cluster_statistic(w,m);

        flip_cluster(w,rng);
        if((i+1)%1000==0) {
            printf("themral : %d\n",i+1);
        }
    }
    printf("end of thermalization!\n");

    // measurement
    double dt = T/100.0;
    int ntime = (int)(T/dt+1);
    double* time_list = (double*)malloc(sizeof(double)*ntime);
    for(int i=0;i<ntime;i++) {
        time_list[i] = (dt*i)/T;
    }

    int ntrial=0;
    for(int i_sweep=0;i_sweep<nsweep;) {
        for(int i=0;i<nskip;i++) {
            remove_vertices(w);
            swapping_graphs(w,m,rng);
            insert_vertices(w,m,rng);
            boundary_condition_initial_state(w,m,initial_condition_type,rng);
            boundary_condition_final_state(w,m,pnif,final_condition_type,rng);
            clustering(w,m);

            if((i+1)==nskip){
                cluster_statistic(w,m);
            }

            flip_cluster(w,rng);
        }
        ntrial++;

        if((ninfected_initial_state(w)==1 && ninfected_final_state(w)>nif) || nocheck_for_measurement) {
            ntrial_ave+=ntrial;
            measurement(w,m,time_list,ntime,block_size);
            i_sweep++;

            ntrial=0;
        }
    }

    // print the sanpshot of final state
    if(0) {
        char snapshot_filename[128];
        for(int i=0;i<10;i++){
            remove_vertices(w);
            swapping_graphs(w,m,rng);
            insert_vertices(w,m,rng);
            boundary_condition_initial_state(w,m,initial_condition_type,rng);
            boundary_condition_final_state(w,m,pnif,final_condition_type,rng);
            clustering(w,m);
            flip_cluster(w,rng);

            sprintf(snapshot_filename,"snapshot_%d.out",i);
            FILE* snapshot_file = fopen(snapshot_filename,"w");

            snapshot_show(w,m,snapshot_file);

            fclose(snapshot_file);
        }
    }

    // free memory
    free(infected_ratio);
    free(time_list);
    free_world_line(w);
    free_model(m);
    gsl_rng_free(rng);
    free(edges);
}
