#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <omp.h>

#include "dtype.h"
#include "union_find.h"

//#define using_omp

static double d_ave=-1;
static double d_max=-1;
static void generate_uniform_sequence(world_line_omp* w, model* m, gsl_rng** rng) {
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

    double lam = (m->sweight)*(w->beta)*d_max/d_ave/(w->nthread);

#ifdef using_omp
    omp_set_num_threads(w->nthread);
    #pragma omp parallel
    {
        int i_thread = omp_get_thread_num();
#else
    for(int i_thread=0;i_thread<(w->nthread);i_thread++) {
#endif
        double* insert_seq = w->insert_seq[i_thread];
        int* insert_bond = w->insert_bond[i_thread];
        int insert_cap = w->insert_cap[i_thread];
        
        double k=0;
        int n=0;
        int bond;
        double dis = gsl_rng_uniform_pos(rng[i_thread]);
        k -= log(dis)/lam;
        while((k<1.0) && (n<insert_cap)) {
            bond = gsl_rng_uniform_pos(rng[i_thread])*(m->nbond);
            if(gsl_rng_uniform_pos(rng[i_thread])*d_max<(m->bond2weight[bond])){
                insert_seq[n]  = k;
                insert_bond[n] = bond;
                n++;
            }

            dis = gsl_rng_uniform_pos(rng[i_thread]);
            k -= log(dis)/lam;
        }
        while(k<1.0) {
            dis = gsl_rng_uniform_pos(rng[i_thread]);
            k -= log(dis)/lam;
            n++;
        }
        w->insert_len[i_thread] = n;
        if(n>insert_cap) {
            w->insert_len[i_thread] = insert_cap;
            realloc_world_line_omp_insert(w,n*2,i_thread);
        }
    }
}

void remove_vertices_omp(world_line_omp* w) {
#ifdef using_omp
    omp_set_num_threads(w->nthread);
    #pragma omp parallel
    {
        int i_thread = omp_get_thread_num();
#else
    for(int i_thread=0;i_thread<(w->nthread);i_thread++) {
#endif
        vertex* v;

        vertex* sequence1 = w->sequenceB[i_thread];
        vertex* sequence2 = w->sequenceA[i_thread];
        if(w->flag[i_thread]) {
            sequence1 = w->sequenceA[i_thread];
            sequence2 = w->sequenceB[i_thread];
        }

        int check_delete;
        int i,j,k;
        k=0;
        for(i=0;i<(w->len[i_thread]);i++) {
            v = &(sequence1[i]);
            check_delete=1;

            for(j=0;j<(v->hNspin);j++) {
                if(v->state[j] != v->state[j+v->hNspin])
                    check_delete=0;
            }

            if(!check_delete) {
                copy_vertex(&(sequence2[k]),&(sequence1[i]));
                k++;
            }
        }
        w->len[i_thread] = k;
        w->flag[i_thread] = !(w->flag[i_thread]);
    }
}

void insert_vertices_omp(world_line_omp* w, model* m, gsl_rng** rng) {
    generate_uniform_sequence(w,m,rng);

    for(int i_thread=0;i_thread<(w->nthread);i_thread++) {
        int len = (w->len[i_thread]) + (w->insert_len[i_thread]);
        realloc_world_line_omp_vertex(w,len,i_thread);        
    }

#ifdef using_omp
    omp_set_num_threads(w->nthread);
    #pragma omp parallel
    {
        int i_thread = omp_get_thread_num();
#else
    for(int i_thread=0;i_thread<(w->nthread);i_thread++) {
#endif

        vertex* sequence1 = w->sequenceB[i_thread];
        vertex* sequence2 = w->sequenceA[i_thread];
        if(w->flag[i_thread]) {
            sequence1 = w->sequenceA[i_thread];
            sequence2 = w->sequenceB[i_thread];
        }

        int nsite = w->nsite;
        int* pstate = &(w->pstate[i_thread*nsite]);

        for(int i=0;i<nsite;i++) 
            pstate[i] = w->istate[i_thread*nsite+i];

        vertex* v;
        int n,i,k,i_site,index;
        double tau1,tau2;

        int mhnspin = m->mhnspin;
        int lstate[mhnspin];

        double* insert_seq = w->insert_seq[i_thread];
        int* insert_bond = w->insert_bond[i_thread];

        k=0;
        n=0;
        tau1=0;
        if(w->len[i_thread]>0) tau1 = (sequence1[0]).tau;

        for(i=0;i<(w->insert_len[i_thread]);i++) {
            tau2 = insert_seq[i];

            while((tau1<tau2) && (k<(w->len[i_thread]))) {
                v = &(sequence1[k]);
                for(i_site=0;i_site<(v->hNspin);i_site++) {
                    index = m->bond2index[v->bond*mhnspin+i_site];
                    pstate[index] = v->state[v->hNspin+i_site];
                }

                copy_vertex(&(sequence2[n]),v);
                n++;
                k++;

                if(k<(w->len[i_thread])) {
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

        while(k<(w->len[i_thread])) {
            copy_vertex(&(sequence2[n]),&(sequence1[k]));
            n++;
            k++;
        }

        w->len[i_thread] = n;
        w->flag[i_thread] = !(w->flag[i_thread]);
    }
}

void clustering_inner_omp(world_line_omp* w, model* m) {

    realloc_world_line_omp_cluster(w);

#ifdef using_omp
    omp_set_num_threads(w->nthread);
    #pragma omp parallel
    {
        int i_thread = omp_get_thread_num();
#else
    for(int i_thread=0;i_thread<(w->nthread);i_thread++) {
#endif
        
        vertex* v;
        int bond,hNspin,t,idn,idp;
        int i,j,index;
        int mnspin = w->mnspin;
        int nsite = w->nsite;
        int* rule;
        int* indices;

        int* first = &(w->first[nsite*i_thread]);
        int* last  = &(w->last[nsite*i_thread]);

        for(i=0;i<nsite;i++) {
            last[i]  = -1;
            first[i] = -1;
        }

        vertex* sequence = w->sequenceB[i_thread];
        if(w->flag) 
            sequence = w->sequenceA[i_thread];

        for(i=0;i<(w->len[i_thread]);i++) {
            v       = &(sequence[i]);
            bond    = v->bond;
            hNspin  = v->hNspin;
            t       = m->bond2type[bond];
            rule    = &(m->link[4*(m->mhnspin)*t]);
            indices = &(m->bond2index[bond*(m->mhnspin)]);
            
            for(j=0;j<2*hNspin;j++) {
                idn = i_thread*(w->mcap)*mnspin+i*mnspin+j;
                idp = i_thread*(w->mcap)*mnspin+i*mnspin+rule[j];
                w->cluster[idn] = idp;
                w->weight[idn]  = rule[2*hNspin+j];
            }

            for(j=0;j<hNspin;j++) {
                index = indices[j];
                idp = i_thread*(w->mcap)*mnspin+i*mnspin+j;
                idn = i_thread*(w->mcap)*mnspin+i*mnspin+j+hNspin;
                if(first[index]==-1) {
                    first[index] = idp;
                    last[index]  = idn;
                } else {
                    merge(w->cluster,w->weight,last[index],idp);
                    last[index] = idn;
                }
            }
        }
    }
}

//without parallelization
void clustering_crossing(world_line_omp* w) {
    int nsite = w->nsite;
    int* first_all = w->first;
    int* last_all = w->last;

    for(int i_thread=1;i_thread<(w->nthread);i_thread++) {
        int* first = &(w->first[nsite*i_thread]);
        int* last  = &(w->last[nsite*i_thread]);

        for(int i=0;i<nsite;i++) {
            if(first_all[i]==-1){
                first_all[i] = first[i];
                last_all[i] = last[i];
            } else if(first[i]==-1) {
                first[i] = last_all[i];
                last[i] = last_all[i];
            } else {
                merge(w->cluster,w->weight,last_all[i],first[i]);
                last_all[i] = last[i];
            }
        }
    }

    for(int i=0;i<nsite;i++) {
        if(first_all[i]!=-1) {
            merge(w->cluster,w->weight,first_all[i],last_all[i]);
        }
    }
}

void flip_cluster_omp(world_line_omp* w, gsl_rng** rng) {
    
#ifdef using_omp
    omp_set_num_threads(w->nthread);
    #pragma omp parallel
    {
        int i_thread = omp_get_thread_num();
#else
    for(int i_thread=0;i_thread<(w->nthread);i_thread++) {
#endif
        
        int start = i_thread*(w->mcap)*(w->mnspin);
        int end = start+(w->len[i_thread])*(w->mnspin);
        double dis;

        for(int i=start;i<end;i++) {
            if((w->cluster[i])==i) {
                dis = gsl_rng_uniform_pos(rng[i_thread]);
                if(dis<0.5) {
                    w->weight[i] = 1;
                } else {
                    w->weight[i] = -1;
                }
            }
        }
#ifdef using_omp
        #pragma omp barrier
#endif

        int* state;
        int hNspin,idv,idr,i,j;

        int mnspin = w->mnspin;

        vertex* sequence = w->sequenceB[i_thread];
        if(w->flag[i_thread])
            sequence = w->sequenceA[i_thread];

        for(i=0;i<(w->len[i_thread]);i++) {
            hNspin = (sequence[i]).hNspin;
            state  = (sequence[i]).state;

            for(j=0;j<2*hNspin;j++) {
                idv = start+i*mnspin+j;
                idr = root(w->cluster,idv);
                
                state[j] = state[j]*(w->weight[idr]);
            }
        }
    }

    for(int i_thread=0;i_thread<(w->nthread);i_thread++) {
        int id,i,j,t,p;
        int mnspin = w->mnspin;
        int mcap = w->mcap;
        int nsite = w->nsite;

        vertex* sequence;

        for(i=0;i<nsite;i++) {
            id = w->first[i+i_thread*nsite];
            if(id!=-1) {
                t  = id/(mnspin*mcap);
                id = id%(mnspin*mcap);
                p  = id/mnspin;
                j  = id%mnspin;

                sequence = w->sequenceB[t];
                if(w->flag[i_thread])
                    sequence = w->sequenceA[t];

                w->istate[i+nsite*i_thread] = (sequence[p]).state[j];
            } else if(i_thread==0) {
                if(gsl_rng_uniform_pos(rng[0])<0.5) {
                    w->istate[i] =  1;
                } else {
                    w->istate[i] = -1;
                }
            } else {
                w->istate[i+nsite*i_thread] = w->istate[i];
            }
        }
    }
}
