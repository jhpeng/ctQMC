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

int ninfected_initial_state(world_line* w) {
    int nnode=w->nsite;
    int ninfected=0;

    for(int i=0;i<nnode;i++) {
        ninfected+=w->istate[i];
    }

    return (ninfected+nnode)/2;
}

int* boundary_condition_frozen_list=NULL;
void boundary_condition_initial_state(world_line* w, model* m, int type, gsl_rng* rng) {
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

    if(boundary_condition_frozen_list==NULL) {
        boundary_condition_frozen_list = (int*)malloc(sizeof(int)*nnode);
    }

    int n=0;
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
        for(int i=0;i<nnode;i++) {
            if(w->istate[i]==-1) {
                boundary_condition_frozen_list[i]=1;
            } else {
                boundary_condition_frozen_list[i]=0;
                ninfected++;
            }
        }
        if(ninfected==1) {
            int check=1;
            while(check) {
                int j = gsl_rng_uniform_pos(rng)*nnode;
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
    int nnode = m->nsite;
    int nbond = m->nbond;
    int length = nnode+(w->nvertices);
    realloc_world_line(w,length);

    vertex* sequence = w->sequenceB;
    if(w->flag) {
        sequence = w->sequenceA;
    }

    int* pstate = w->pstate;
    int n=w->nvertices;

    if(type==0) {
        for(int i=0;i<nnode;i++) {
            (sequence[n]).tau      = 1.0;
            (sequence[n]).bond     = nbond+i;
            (sequence[n]).hNspin   = 1;
            (sequence[n]).state[0] = pstate[i];
            (sequence[n]).state[1] = pstate[i];
            n++;
        }
    } else if(type==1) {
        int inf=0;
        for(int i=0;i<nnode;i++) inf += (pstate[i]+1)/2;

        double pdis = 1.0;
        if(inf!=0) pdis=(p*nnode)/inf;
        
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


    w->nvertices = n;
}

static void print_state(int* state, int nnode) {
    for(int i=0;i<nnode;i++) {
        printf("%d ",(state[i]+1)/2);
    }
    printf("\n");
}

static void save_state(FILE* file, int* state, int nnode) {
    for(int i=0;i<nnode;i++) {
        fprintf(file,"%d ",(state[i]+1)/2);
    }
    fprintf(file,"\n");
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

time_t start_time,end_time;
unsigned long int measurement_count=0;
double* infected_ratio=NULL;
double* infected_time=NULL;
double total_infected_time=0;
void measurement(world_line* w, model* m, double* time_list, int ntime, int block_size) {
    if(infected_ratio==NULL) {
        infected_time = (double*)malloc(sizeof(double)*(w->nsite));
        infected_ratio = (double*)malloc(sizeof(double)*ntime);
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

        double ninfection = ninfection_ave_value();
        double nrecover  = nrecover_ave_value();
        total_infected_time = total_infected_time/block_size*(w->beta);
        fprintf(file_g,"%.12e %.12e %.12e\n",ninfection,nrecover,total_infected_time);

        print_ncluster_flippable();
        print_ninfection();
        print_nrecover();
        printf("total infected time = %.12e\n",total_infected_time);

        save_configuration(file_conf,w,m,time_list,ntime);

        total_infected_time=0;
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
    int block_size=atoi(argv[5]);
    int nblock = atoi(argv[6]);
    int thermal = atoi(argv[7]);
    int nsweep  = nblock*block_size;
    unsigned long int seed=atoi(argv[8]);

    int nnode;
    int nedge;
    int* edges = read_edgelist(filename,&nnode,&nedge);
    double pnif = ((double)nif)/nnode;


    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);


    model* m = sis_model_uniform_infection(alpha,gamma,nnode,nedge,edges);
    world_line* w = malloc_world_line(1024,2*(m->mhnspin),m->nsite);

    // initial state
    for(int i=0;i<(w->nsite);i++) {
        w->istate[i] = -1;
        //if(gsl_rng_uniform_pos(rng)<0.05) {
        if(i<nif) {
            w->istate[i] = 1;
        }
    }

    // thermalization
    w->beta = T;
    for(int i=0;i<thermal;i++) {
        remove_vertices(w);
        swapping_graphs(w,m,rng);
        insert_vertices(w,m,rng);
        boundary_condition_initial_state(w,m,0,rng);
        boundary_condition_final_state(w,m,pnif,2,rng);
        clustering(w,m);
        flip_cluster(w,rng);
        if((i+1)%block_size==0) {
            printf("themral : %d\n",i+1);
        }
    }
    printf("end of thermalization!\n");

    // measurement
    double dt = 0.1;
    int ntime = (int)(T/dt+1);
    double* time_list = (double*)malloc(sizeof(double)*ntime);
    for(int i=0;i<ntime;i++) {
        time_list[i] = (dt*i)/T;
    }

    for(int i_sweep=0;i_sweep<nsweep;) {
        for(int i=0;i<10;i++) {
            remove_vertices(w);
            swapping_graphs(w,m,rng);
            insert_vertices(w,m,rng);
            boundary_condition_initial_state(w,m,0,rng);
            boundary_condition_final_state(w,m,pnif,2,rng);
            clustering(w,m);
            flip_cluster(w,rng);
        }

        //if(ninfected_initial_state(w)==1) {
        if(1) {
            measurement(w,m,time_list,ntime,block_size);
            i_sweep++;
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
