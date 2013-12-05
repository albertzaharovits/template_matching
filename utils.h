#ifndef UTILS_H
#define UTILS_H

#include "colorimage.h"
#include "settings.h"

#define LMAGENTA 60.31993f
#define AMAGENTA 98.25421f
#define BMAGENTA -60.84298f

namespace Utils {

  static const fp PI = 3.14159265;
  static const fp D2R = PI/180.0;

  class CircularSamplingData {
      public:
        uint cis_start;
        uint cis_step;
        uint cis_n;
        float* cis_l;
        float* cis_a;
        float* cis_b;

        CircularSamplingData(uint cis_start, uint cis_step, uint cis_n);
        ~CircularSamplingData();
        CircularSamplingData(CircularSamplingData&& other);
        CircularSamplingData& operator=(CircularSamplingData&& other);
        CircularSamplingData(const CircularSamplingData& other);
  };

  void circle_pix_mean( unsigned int yc, unsigned int xc, unsigned int r,
                        const Image::ColorImage& im, fp& _l, fp& _a, fp& _b);


  CircularSamplingData circle_sampling( const Image::ColorImage& im,
                                        uint circle_start, uint circle_step_delta);
}

#endif // UTILS_H
