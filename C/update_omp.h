#ifndef update_omp_h
#define update_omp_h

#include <gsl/gsl_rng.h>

#include "dtype.h"

void remove_vertices_omp(world_line_omp* w);

void insert_vertices_omp(world_line_omp* w, model* m, gsl_rng** rng);

void clustering_inner_omp(world_line_omp* w, model* m);

void clustering_crossing(world_line_omp* w);

void clustering_crossing_omp(world_line_omp* w);

void flip_cluster_omp(world_line_omp* w, gsl_rng** rng);

int check_world_line_omp_configuration(world_line_omp* w, model* m);

#endif
