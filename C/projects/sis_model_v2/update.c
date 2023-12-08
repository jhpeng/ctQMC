#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "union_find.h"

//#define weighted_sampling

#ifdef weighted_sampling
static int weighted_sampling(double* cmf, int length, gsl_rng* rng) {
    double dis = gsl_rng_uniform_pos(rng)*cmf[length-1];

    int i = -1;
    int j = length;
    int d = j-1;

    while(d>1) {
        if(dis>cmf[i+d/2]){
            i = i+d/2;
        } else {
            j = i+d/2;
        }
        d = j-i;
    }

    return j;
}
#endif

static double* insert_seq;
static int* insert_bond;
static int insert_len=0;
static int insert_cap=0;
static double d_ave=-1;
static double d_max=-1;
static void uniform_sequence_sampling(model* m, double lam, double start, gsl_rng* rng) {
    if(d_max<0){
        d_max = 0;
        d_ave = 0;
        for(int i=0;i<(m->nbond);i++) {
            d_ave += m->bond2weight[i];
            if((m->bond2weight[i])>d_max)
                d_max = m->bond2weight[i];
        }
        d_ave = d_ave/(m->nbond);
    }

    if(insert_cap==0) {
        insert_cap = (int)(lam+sqrt(lam)*10+1024);
        insert_seq  = (double*)malloc(sizeof(double)*insert_cap);
        insert_bond = (int*)malloc(sizeof(int)*insert_cap);
    }

    lam = lam*d_max/d_ave;

    double k=0;
    int n=0;
    int bond;

    double dis = gsl_rng_uniform_pos(rng);
    k -= log(dis)/lam;
    while((k<1.0) && (n<insert_cap)) {
        bond = gsl_rng_uniform_pos(rng)*(m->nbond);
        if(gsl_rng_uniform_pos(rng)*d_max<(m->bond2weight[bond])){
            insert_seq[n]  = k+start;
            insert_bond[n] = bond;
            n++;
        }

        dis = gsl_rng_uniform_pos(rng);
        k -= log(dis)/lam;
    }
    while(k<1.0) {
        dis = gsl_rng_uniform_pos(rng);
        k -= log(dis)/lam;
        n++;
    }
    insert_len = n;
    if(n>insert_cap) {
        double* seq = (double*)malloc(sizeof(double)*n*2);
        int* bond   = (int*)malloc(sizeof(int)*n*2);
        for(int i=0;i<insert_cap;i++){
            seq[i]  = insert_seq[i];
            bond[i] = insert_bond[i];
        }
        free(insert_seq);
        insert_seq  = seq;
        insert_bond = bond;
        insert_len = insert_cap;
        insert_cap = n*2;
    }
}

static int ninfection=0;
static int nrecover=0;

double ninfection_value() {
    return ninfection;
}

double nrecover_value() {
    return nrecover;
}

void remove_vertices(world_line* w) {
    vertex* v;

    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    ninfection=0;
    nrecover=0;

    int check_delete;
    int i,j,k;
    k=0;
    for(i=0;i<w->nvertices;i++) {
        v = &(sequence1[i]);
        check_delete = 1;

        for(j=0;j<(v->hNspin);j++) {
            if(v->state[j]!=v->state[j+v->hNspin])
                check_delete = 0;
        }

        if(!check_delete) {
            copy_vertex(&(sequence2[k]),&(sequence1[i]));
            k++;

            if(v->hNspin==1) {
                nrecover++;
            } else if(v->hNspin==2) {
                ninfection++;
            }
        }
    }
    w->nvertices = k;
    w->flag = !(w->flag);
}

void swapping_graphs(world_line* w, model* m, gsl_rng* rng) {
    int nnode = m->nsite;
    int nedge = (m->nbond-3*nnode)/7;
    vertex* v;

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    for(int i=0;i<(w->nvertices);i++) {
        v = &(sequence[i]);

        // swap between type (1,3,5) or (2,4,6)
        int type = m->bond2type[v->bond];
        if((type==1 || type==3) || type==5) {
            int i_edge = (v->bond) % nedge;
            double random_value = gsl_rng_uniform_pos(rng)*3.0;
            if(random_value<1.0) {
                v->bond = i_edge + nedge*1;
            } else if(random_value<2.0) {
                v->bond = i_edge + nedge*3;
            } else {
                v->bond = i_edge + nedge*5;
            }
        } else if((type==2 || type==4) || type==6) {
            int i_edge = (v->bond) % nedge;
            double random_value = gsl_rng_uniform_pos(rng)*3.0;
            if(random_value<1.0) {
                v->bond = i_edge + nedge*2;
            } else if(random_value<2.0) {
                v->bond = i_edge + nedge*4;
            } else {
                v->bond = i_edge + nedge*6;
            }
        }

        // swap between type (7,8)
        else if(type==7) {
            if(gsl_rng_uniform_pos(rng)<0.5) {
                v->bond += nnode;
            }
        } else if(type==8) {
            if(gsl_rng_uniform_pos(rng)<0.5) {
                v->bond -= nnode;
            }
        }
    }

}

void insert_vertices(world_line* w, model* m, gsl_rng* rng) {
    double lam = (m->sweight)*(w->beta);

    uniform_sequence_sampling(m,lam,0,rng);

    int length = insert_cap+w->nvertices;
    realloc_world_line(w,length);

    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    int* pstate = w->pstate;
    int  nsite = w->nsite;

    for(int i=0;i<nsite;i++) pstate[i] = w->istate[i];

    vertex* v;
    int n,i,k,i_site,index;
    double tau1,tau2;

    int mhnspin = m->mhnspin;
    int lstate[mhnspin];

    k=0;
    n=0;

    tau1 = 0;
    if(w->nvertices!=0) tau1 = (sequence1[0]).tau;

    for(i=0;i<insert_len;i++) {
        tau2 = insert_seq[i];

        while((tau1<tau2) && (k<(w->nvertices))) {
            v = &(sequence1[k]);
            for(i_site=0;i_site<(v->hNspin);i_site++) {
                index = m->bond2index[v->bond*mhnspin+i_site];
#if 0
                if(index<0 || index>=nsite) {
                    printf("index = %d \n",index);
                    printf("------- vertex information ------- \n");
                    printf("tau    = %.16lf \n",v->tau);
                    printf("bond   = %d \n",v->bond);
                    printf("hNspin (w) = %d \n",v->hNspin);
                    printf("hNspin (m) = %d \n",m->bond2hNspin[v->bond]);

                    printf("state (");
                    for(int i_state=0;i_state<2*(v->hNspin);i_state++) {
                        printf(" %d",v->state[i_state]);
                    }
                    printf(")\n");

                    printf("indices (");
                    for(int i_state=0;i_state<(v->hNspin);i_state++) {
                        printf(" %d",m->bond2index[v->bond*mhnspin+i_state]);
                    }
                    printf(")\n");

                    exit(-1);
                }
#endif
                pstate[index] = v->state[v->hNspin+i_site];
            }

            copy_vertex(&(sequence2[n]),v);
            n++;
            k++;

            if(k<(w->nvertices)){
                tau1 = (sequence1[k]).tau;
            }
        }

        if(tau1!=tau2) {
            int bond     = insert_bond[i];
            int t        = m->bond2type[bond];
            int hNspin   = m->bond2hNspin[bond];
            insert_rule rule = m->insert[t];

            for(i_site=0;i_site<hNspin;i_site++) { 
                index = m->bond2index[bond*mhnspin+i_site];
                lstate[i_site] = pstate[index];
            }

            if(rule(lstate)) {
                (sequence2[n]).tau    = tau2;
                (sequence2[n]).bond   = bond;
                (sequence2[n]).hNspin = hNspin;

                for(i_site=0;i_site<hNspin;i_site++) {
                    (sequence2[n]).state[i_site]        = lstate[i_site];
                    (sequence2[n]).state[i_site+hNspin] = lstate[i_site];
                }

                n++;
            }
        }
    }

    while(k<(w->nvertices)) {
        copy_vertex(&(sequence2[n]),&(sequence1[k]));
        n++;
        k++;
    }

    w->nvertices = n;
    w->flag = !(w->flag);
}

void clustering(world_line* w, model* m) {
    vertex* v;
    int bond,hNspin,t,idn,idp;
    int i,j,index;
    int mnspin = w->mnspin;
    int nsite = w->nsite;
    int* rule;
    int* indices;

    int* first = w->first;
    int* last  = w->last;

    for(i=0;i<nsite;i++) {
        last[i]  = -1;
        first[i] = -1;
    }

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    for(i=0;i<(w->nvertices);i++) {
        v       = &(sequence[i]);
        bond    = v->bond;
        hNspin  = v->hNspin;
        t       = m->bond2type[bond];
        rule    = &(m->link[4*(m->mhnspin)*t]);
        indices = &(m->bond2index[bond*(m->mhnspin)]);

        for(j=0;j<2*hNspin;j++) {
            idn = i*mnspin+j;
            idp = i*mnspin+rule[j];
            w->cluster[idn] = idp;
            w->weight[idn]  = rule[2*hNspin+j];
        }

        for(j=0;j<hNspin;j++) {
            index = indices[j];
            idp = i*mnspin+j;
            idn = i*mnspin+j+hNspin;
            if(first[index]==-1) {
                first[index] = idp;
                last[index]  = idn;
            } else {
                merge(w->cluster,w->weight,last[index],idp);
                last[index] = idn;
            }
        }
    }

/*  disable for open boundary
**  for(i=0;i<nsite;i++) {
**      if(first[i]!=-1) {
**          merge(w->cluster,w->weight,first[i],last[i]);
**      }
**  }
*/
}

int  cluster_statistic_length=0;
int  cluster_statistic_counter=0;
int* cluster_statistic_count=NULL;
int* cluster_statistic_fcluster=NULL;
int* cluster_statistic_infection=NULL;
double* cluster_statistic_taus=NULL;
double* cluster_statistic_infection_size=NULL;
void cluster_statistic(world_line* w, model* m) {
    int idv,idr,i,j,index;
    int bond,hNspin;
    double tau;
    int* indices;
    int* state;
    vertex* v = NULL;

    int nsite = w->nsite;
    int mnspin = w->mnspin;
    int length = w->length;

    if(cluster_statistic_taus==NULL) {
        cluster_statistic_taus = (double*)malloc(sizeof(double)*nsite);

        if(cluster_statistic_taus==NULL) {
            printf("Memory Allocating Error : update.c (cluster_statistic)\n");
            exit(-1);
        }

        cluster_statistic_infection_size = (double*)malloc(sizeof(double)*nsite);

        if(cluster_statistic_infection_size==NULL) {
            printf("Memory Allocating Error : update.c (cluster_statistic)\n");
            exit(-1);
        }

        cluster_statistic_fcluster = (int*)malloc(sizeof(int)*nsite);

        if(cluster_statistic_fcluster==NULL) {
            printf("Memory Allocating Error : update.c (cluster_statistic)\n");
            exit(-1);
        }

        cluster_statistic_infection = (int*)malloc(sizeof(int)*nsite);

        if(cluster_statistic_infection==NULL) {
            printf("Memory Allocating Error : update.c (cluster_statistic)\n");
            exit(-1);
        }
    }

    int* fcluster = cluster_statistic_fcluster;
    int* isize = cluster_statistic_infection;
    for(i=0;i<nsite;i++) {
        fcluster[i] = 0;
        isize[i] = 0;
    }

    if(length*mnspin>cluster_statistic_length) {
        cluster_statistic_length = length*mnspin;
        if(cluster_statistic_count!=NULL)
            free(cluster_statistic_count);

        cluster_statistic_count = (int*)malloc(sizeof(int)*cluster_statistic_length);

        if(cluster_statistic_count==NULL) {
            printf("Memory Allocating Error : update.c (cluster_statistic)\n");
            exit(-1);
        }
    }

    for(i=0;i<cluster_statistic_length;i++) {
        cluster_statistic_count[i] = 1;
    }

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    int number_of_free_cluster=0;
    int number_of_cluster=0;
    int size_of_free_cluster=0;
    int size_of_cluster=0;
    double cluster_size_in_time=0;
    double infection_size_in_time=0;
    for(i=0;i<(w->nvertices);i++) {
        v       = &(sequence[i]);
        tau     = v->tau;
        bond    = v->bond;
        hNspin  = v->hNspin;
        indices = &(m->bond2index[bond*(m->mhnspin)]);
        state   = v->state;

        for(j=0;j<2*hNspin;j++) {
            idv = i*mnspin+j;
            idr = root(w->cluster,idv);
            if(cluster_statistic_count[idr]) {
                number_of_cluster++;

                if(w->weight[idr]>0) {
                    number_of_free_cluster++;
                    size_of_free_cluster += w->weight[idr];
                    size_of_cluster += w->weight[idr];
                } else {
                    size_of_cluster -= w->weight[idr];
                }

                cluster_statistic_count[idr]=0;
            }
        }

        for(j=0;j<hNspin;j++) {
            index = indices[j];
            idv = i*mnspin+j+hNspin;

            if(fcluster[index]) {
                cluster_size_in_time += (tau - cluster_statistic_taus[index]);
                fcluster[index] = 0;
            }

            idr = root(w->cluster,idv);
            if(w->weight[idr]>0) {
                cluster_statistic_taus[index] = tau;
                fcluster[index] = 1;
            }

            if(isize[index]) {
                infection_size_in_time += (tau - cluster_statistic_infection_size[index]);
                isize[index] = 0;
            }

            if(state[j+hNspin]==1) {
                cluster_statistic_infection_size[index] = tau;
                isize[index] = 1;
            }
        }
    }
    
    double ratio1 = ((double)number_of_free_cluster)/((double)number_of_cluster);
    double ratio2 = ((double)size_of_free_cluster)/((double)size_of_cluster);
    cluster_statistic_counter++;
    printf("---------------------------------------------------\n");
    printf("n = %d \n", cluster_statistic_counter);
    printf("number of cluster = %d, number of free cluster = %d, ratio = %.16lf \n",number_of_cluster,number_of_free_cluster,ratio1);
    printf("size of cluster = %d, size of free cluster = %d, ratio = %.16lf \n",size_of_cluster,size_of_free_cluster,ratio2);
    printf("size of cluster (t) = %lf\n",cluster_size_in_time);
    printf("infection size in time (t) = %lf\n",infection_size_in_time);

    char filename[128] = "cluster_statistic.txt";
    FILE* sfile = fopen(filename,"a");
    fprintf(sfile,"%.12e %.12e %d \n", cluster_size_in_time, infection_size_in_time, w->nvertices);
    fclose(sfile);
}

void flip_cluster(world_line* w, gsl_rng* rng) {
    int* state;
    int hNspin,idv,idr,id,p,i,j;

    int mnspin = w->mnspin;
    int nsite  = w->nsite;

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    for(i=0;i<(w->nvertices);i++) {
        hNspin = (sequence[i]).hNspin;
        state  = (sequence[i]).state;

        for(j=0;j<2*hNspin;j++) {
            idv = i*mnspin+j;
            idr = root(w->cluster,idv);
            if(w->weight[idr]>0) {
                if(gsl_rng_uniform_pos(rng)<1.0) {
                    w->weight[idr] =  0;
                } else {
                    w->weight[idr] = -1;
                }
            }
            if(w->weight[idr]==0) {
                state[j] = -state[j];
            }
        }
    }

    for(i=0;i<nsite;i++) {
        id = w->first[i];
        if(id!=-1) {
            p = id/mnspin;
            j  =id%mnspin;
            w->istate[i] = (sequence[p]).state[j];
        } else if(gsl_rng_uniform_pos(rng)<0.5) {
            w->istate[i] =  1;
        } else {
            w->istate[i] = -1;
        }

        id = w->last[i];
        if(id!=-1) {
            p = id/mnspin;
            j  =id%mnspin;
            w->pstate[i] = (sequence[p]).state[j];
        }
    }
}

void remove_only_fixed_vertices(world_line* w) {
    int hNspin,idv,idr,i,j,k;

    int mnspin = w->mnspin;

    vertex* v;

    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    ninfection=0;
    nrecover=0;

    int check_save;
    k=0;
    for(i=0;i<w->nvertices;i++) {
        v = &(sequence1[i]);
        hNspin = v->hNspin;
        check_save = 0;

        for(j=0;j<hNspin;j++) {
            if(v->state[j]!=v->state[j+v->hNspin]) {
                check_save = 1;
            } 
        }

        for(j=0;j<2*hNspin;j++) {
            idv = i*mnspin+j;
            idr = root(w->cluster,idv);
            if(w->weight[idr]==0) {
                check_save = 1;
            }
        }

        if(check_save) {
            copy_vertex(&(sequence2[k]),&(sequence1[i]));
            k++;

            if(v->hNspin==1) {
                nrecover++;
            } else if(v->hNspin==2) {
                ninfection++;
            }
        }
    }
    w->nvertices = k;
    w->flag = !(w->flag);
}

void snapshot_show(world_line* w, model* m, FILE* file) {
    vertex* v;
    int idv,idr;
    int bond,hNspin,type,index;
    int* indices;
    double tau;

    int nsite = w->nsite;
    int mnspin = w->mnspin;

    vertex* sequence = w->sequenceB;
    if(w->flag) 
        sequence = w->sequenceA;

    printf("------------------ snapshot -------------------\n");
    printf("data format: (type tau [node_id, ] free_cluster) \n");

    for(int i=0;i<nsite;i++) {
        if(w->istate[i]>0) {
            type  = 10;
            tau   = 0.0;
            index = i;
            printf("%d %f %d 0 \n",type,tau,index);
            fprintf(file,"%d %f %d 0 \n",type,tau,index);
        }
    }

    for(int i=0;i<(w->nvertices);i++) {
        v       = &(sequence[i]);
        bond    = v->bond;
        hNspin  = v->hNspin;
        tau     = v->tau;
        type    = m->bond2type[bond];
        indices = &(m->bond2index[bond*(m->mhnspin)]);

        int check_condition = 0;

        for(int j=0;j<hNspin;j++) {
            if(v->state[j]!=v->state[j+v->hNspin]) {
                check_condition = 1;
            } 
        }

        for(int j=0;j<2*hNspin;j++) {
            idv = i*mnspin+j;
            idr = root(w->cluster,idv);
            if(w->weight[idr]==0) {
                check_condition = 1;
            }
        }

        if(check_condition) {
            printf("%d %f ",type,tau);
            fprintf(file,"%d %f ",type,tau);

            for(int j=0;j<hNspin;j++) {
                index = indices[j];
                printf("%d ",index);
                fprintf(file,"%d ",index);
            }

            int free_cluster=0;
            for(int j=0;j<2*hNspin;j++) {
                idv = i*mnspin+j;
                idr = root(w->cluster,idv);
                if(w->weight[idr]==0) {
                    free_cluster=1;
                }
            }
            printf("%d \n",free_cluster);
            fprintf(file,"%d \n",free_cluster);
        }
    }

    for(int i=0;i<nsite;i++) {
        int free_cluster=0;
        idv = w->last[i];
        if(idv!=-1) {
            idr = root(w->cluster,idv);
            if(w->weight[idr]==0) {
                free_cluster=1;
            }
        }

        if(w->pstate[i]>0) {
            type = 11;
            tau  = 1.0;
            index = i;

            printf("%d %f %d %d \n",type,tau,index,free_cluster);
            fprintf(file,"%d %f %d %d \n",type,tau,index,free_cluster);
        } else if(free_cluster) {
            type = 11;
            tau  = 1.0;
            index = i;

            printf("%d %f %d %d \n",type,tau,index,free_cluster);
            fprintf(file,"%d %f %d %d \n",type,tau,index,free_cluster);
        }
    }
}

