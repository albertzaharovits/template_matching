#include "utils.h"

void Utils::circle_pix_mean( unsigned int yc, unsigned int xc, unsigned int r,
                             const Image::ColorImage& im, fp& _l, fp& _a, fp& _b) {

  unsigned int count = 0;
  _l = _a = _b = 0.0;
  float p;
  unsigned int x, y;

  /* return 0 if no complete circle fits the image */
  if( (xc < r) || ((xc+r) >= im.get_width()) || (yc < r) || ((yc+r) >= im.get_height()) )
    return;

  if( r == 0) {
    count = 1;
    _l = static_cast<fp>( im.L(yc, xc));
    _a = static_cast<fp>( im.A(yc, xc));
    _b = static_cast<fp>( im.B(yc, xc));
  }
  else {
    x = 0;
    y = r;
    /* cross-tip points */
    _l = static_cast<fp>( im.L( r+yc, xc) + im.L( -r+yc, xc) + im.L( yc, r+xc) + im.L( yc, -r+xc) );
    _a = static_cast<fp>( im.A( r+yc, xc) + im.A( -r+yc, xc) + im.A( yc, r+xc) + im.A( yc, -r+xc) );
    _b = static_cast<fp>( im.B( r+yc, xc) + im.B( -r+yc, xc) + im.B( yc, r+xc) + im.B( yc, -r+xc) );
    count = 4;

    p = 1.25f - static_cast<float>( r);

    while(1) {

      p = p + 2.f*static_cast<float>(++x) + 1.f;
      if( p > 0.f)
        p -= 2.f*static_cast<float>(--y);

      if( x >= y) {
        if( x == y) {
          _l += static_cast<fp>( im.L(-y+yc,-x+xc) + im.L(-y+yc, x+xc) +
                                 im.L( y+yc, x+xc) + im.L( y+yc,-x+xc));
          _a += static_cast<fp>( im.A(-y+yc,-x+xc) + im.A(-y+yc, x+xc) +
                                 im.A( y+yc, x+xc) + im.A( y+yc,-x+xc));
          _b += static_cast<fp>( im.B(-y+yc,-x+xc) + im.B(-y+yc, x+xc) +
                                 im.B( y+yc, x+xc) + im.B( y+yc,-x+xc));
          count += 4;
        }
        break;
      }

      /*symmetry points in the other seven octants*/
      _l += static_cast<fp>( im.L( y+yc, x+xc) + im.L( x+yc, y+xc) + im.L(-x+yc, y+xc) + im.L(-y+yc, x+xc) +
                             im.L(-y+yc,-x+xc) + im.L(-x+yc,-y+xc) + im.L( x+yc,-y+xc) + im.L( y+yc,-x+xc) );
      _a += static_cast<fp>( im.A( y+yc, x+xc) + im.A( x+yc, y+xc) + im.A(-x+yc, y+xc) + im.A(-y+yc, x+xc) +
                             im.A(-y+yc,-x+xc) + im.A(-x+yc,-y+xc) + im.A( x+yc,-y+xc) + im.A( y+yc,-x+xc) );
      _b += static_cast<fp>( im.B( y+yc, x+xc) + im.B( x+yc, y+xc) + im.B(-x+yc, y+xc) + im.B(-y+yc, x+xc) +
                             im.B(-y+yc,-x+xc) + im.B(-x+yc,-y+xc) + im.B( x+yc,-y+xc) + im.B( y+yc,-x+xc) );
      count += 8;
    }
  }

  _l = _l/count;
  _a = _a/count;
  _b = _b/count;
}

Utils::CircularSamplingData::CircularSamplingData(uint circle_start, uint circle_step_delta, uint cis_n) {
  this->cis_start = circle_start;
  this->cis_step = circle_step_delta;
  this->cis_n = cis_n;
  posix_memalign( (void**)&cis_l, MEMALLIGN, cis_n*sizeof(float));
  posix_memalign( (void**)&cis_a, MEMALLIGN, cis_n*sizeof(float));
  posix_memalign( (void**)&cis_b, MEMALLIGN, cis_n*sizeof(float));
}

Utils::CircularSamplingData::~CircularSamplingData() {
  this->cis_start = 0;
  this->cis_step = 0;
  this->cis_n = 0;
  if(this->cis_l) {
    free(this->cis_l);
    this->cis_l = NULL;
  }
  if(this->cis_a) {
    free(this->cis_a);
    this->cis_a = NULL;
  }
  if(this->cis_b) {
    free(this->cis_b);
    this->cis_b = NULL;
  }
}

Utils::CircularSamplingData::CircularSamplingData(Utils::CircularSamplingData&& other)
  : cis_start(other.cis_start)
  , cis_step(other.cis_step)
  , cis_n(other.cis_n)
  , cis_l(other.cis_l)
  , cis_a(other.cis_a)
  , cis_b(other.cis_b) {

    other.cis_start = other.cis_step = other.cis_n = 0;
    other.cis_l = other.cis_a = other.cis_b = NULL;
}

Utils::CircularSamplingData::CircularSamplingData(const Utils::CircularSamplingData& other)
  : cis_start(other.cis_start)
  , cis_step(other.cis_step)
  , cis_n(other.cis_n) {

    if(other.cis_l != NULL) {
      posix_memalign( (void**)&cis_l, MEMALLIGN, cis_n*sizeof(float));
      std::copy( other.cis_l, other.cis_l+cis_n, cis_l);
    }
    else
      cis_l = NULL;

    if(other.cis_a != NULL) {
      posix_memalign( (void**)&cis_a, MEMALLIGN, cis_n*sizeof(float));
      std::copy( other.cis_a, other.cis_a+cis_n, cis_a);
    }
    else
      cis_a = NULL;

    if(other.cis_b != NULL) {
      posix_memalign( (void**)&cis_b, MEMALLIGN, cis_n*sizeof(float));
      std::copy( other.cis_b, other.cis_b+cis_n, cis_b);
    }
    else
      cis_a = NULL;
}

Utils::CircularSamplingData& Utils::CircularSamplingData::operator=(Utils::CircularSamplingData&& other) {
  if( (this!=&other) ) {
    if( cis_n!=other.cis_n) {
      cis_start = other.cis_start;
      cis_step = other.cis_step;
      cis_n = other.cis_n;
      if( cis_l != NULL)
        free(cis_l);
      if( cis_a != NULL)
        free(cis_a);
      if( cis_b != NULL)
        free(cis_b);

      posix_memalign( (void**)&cis_l, MEMALLIGN, cis_n*sizeof(float));
      posix_memalign( (void**)&cis_a, MEMALLIGN, cis_n*sizeof(float));
      posix_memalign( (void**)&cis_b, MEMALLIGN, cis_n*sizeof(float));
    }
    std::copy( other.cis_l, other.cis_l+other.cis_n, cis_l);
    std::copy( other.cis_a, other.cis_a+other.cis_n, cis_a);
    std::copy( other.cis_b, other.cis_b+other.cis_n, cis_b);

    other.cis_start = other.cis_step = other.cis_n = 0;
    other.cis_l = other.cis_a = other.cis_b = NULL;
  }

  return *this;
}

Utils::CircularSamplingData Utils::circle_sampling( const Image::ColorImage& im,
                                                    uint circle_start, uint circle_step_delta) {

  uint radius = im.get_radius();
  uint count = (radius-circle_start)/circle_step_delta + 1;
  Utils::CircularSamplingData sdata( circle_start, circle_step_delta, count);

  uint r = circle_start;
  fp l, a, b;

  unsigned int yc = im.get_height()/2;
  unsigned int xc = im.get_width()/2;

  for(uint i=0;i<count;i++) {

    Utils::circle_pix_mean( yc, xc, r, im, l, a, b);
    sdata.cis_l[i] = l;
    sdata.cis_a[i] = a;
    sdata.cis_b[i] = b;
    r += circle_step_delta;
  }

  return sdata;
}
