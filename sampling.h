#ifndef SAMPLING_H
#define SAMPLING_H

#include "settings.h"
#include <stdlib.h>

namespace Sampling {
  class CircularSamplingData {
      public:
        uint id;
        float scale;

        uint cis_start;
        uint cis_step;
        uint cis_n;
        float* cis_l;
        float* cis_a;
        float* cis_b;

        float cis_l_S;
        float cis_l_S2;

        CircularSamplingData( uint cis_start, uint cis_step, uint cis_n, uint id=0, float scale=-33.f);
        ~CircularSamplingData();
        CircularSamplingData(CircularSamplingData&& other);
        CircularSamplingData& operator=(CircularSamplingData&& other);
        CircularSamplingData(const CircularSamplingData& other);

        bool operator<(const CircularSamplingData& other) const;
  };
}

inline
bool Sampling::CircularSamplingData::operator<(const Sampling::CircularSamplingData& other) const {
  return cis_n < other.cis_n;
}

#endif // SAMPLING_H
