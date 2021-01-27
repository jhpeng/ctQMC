#ifndef update_h
#define update_h

void remove_vertices(world_line* w);

void insert_vertices(world_line* w, model* m, gsl_rng* rng);

void clustering(world_line* w, model* m);

void flip_cluster(world_line* w, gsl_rng* rng);

#endif
