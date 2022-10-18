#include <gsl/gsl_rng.h>

#ifndef networks_h
#define networks_h

int* read_edgelist(char* filename, int* nnode, int* nedge);

void nearest_nb_show();

int nearest_nb_random_assign(int i, gsl_rng* rng);

int nearest_nb_arg_max_degree();

#endif
