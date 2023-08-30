#ifndef update_h
#define update_h
#include <gsl/gsl_rng.h>

#include "dtype.h"

double ninfection_value();

double nrecover_value();

void remove_vertices(world_line* w);

void remove_only_fixed_vertices(world_line* w);

void swapping_graphs(world_line* w, model* m, gsl_rng* rng);

void insert_vertices(world_line* w, model* m, gsl_rng* rng);

void clustering(world_line* w, model* m);

void cluster_statistic(world_line* w, model* m);

void flip_cluster(world_line* w, gsl_rng* rng);

//int check_periodic(world_line* w, model* m);

void snapshot_show(world_line* w, model* m, FILE* file);

#endif
