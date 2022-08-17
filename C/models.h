#ifndef models_h
#define models_h

#include "dtype.h"

model* jq3_ladder_square(int lx, int ly, double q);

model* jq3_ladder_square_impurity_spin_half(int lx, int ly, double q);

model* jq3_ladder_square_impurity_spin_one(int lx, int ly, double q);

model* quantum_link_model_2d_square(int lx, int ly, double lambda);

model* quantum_link_model_3d_cubic(int lx, int ly, int lz, double lambda);

#endif
