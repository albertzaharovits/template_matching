#include <algorithm>

#include "utils.h"
#include "colorimage.h"
#include "ControlDict.h"

using namespace std;

/* ColorImage class*/
Image::ColorImage::ColorImage(const std::string file_name) {

  ifstream file(file_name.c_str(), ios::in|ios::binary|ios::ate);
  if( file.is_open()) {

    HeaderStr header;
    file.seekg( 0, ios::beg);
    file.read( (char*)&header, sizeof(HeaderStr));

    set_height( header.height);
    set_width( header.width);
    id = atoi( file_name.substr(0, 3).c_str());

    file.seekg( header.offset, ios::beg);
    int bytes_per_pixel = header.bit_per_pixel / 8;
    int padding = ((header.width*bytes_per_pixel) % 4 == 0) ? 0 : 4 - ((header.width*bytes_per_pixel) % 4);
    /* allocate the size needed to read the BMP data */
    int _size = sizeof(char)*(header.height*(header.width*bytes_per_pixel+padding));
    unsigned char* _data = (unsigned char*) malloc(_size);
    /* read the data */
    file.read((char*)_data, _size);
    posix_memalign( (void**)&l, MEMALLIGN, height*_width_*sizeof(float));
    posix_memalign( (void**)&a, MEMALLIGN, height*_width_*sizeof(float));
    posix_memalign( (void**)&b, MEMALLIGN, height*_width_*sizeof(float));

    unsigned int offset = 0;
    /* in the Bitmap format, pixels are in a reversed order */
    for(int i=header.height-1; i>=0; i--) {
        for(int j=0; j<static_cast<int>(header.width); j++) {

          float var_B = static_cast<float>(_data[offset]) / 255.f;
          float var_G = static_cast<float>(_data[offset+1]) / 255.f;
          float var_R = static_cast<float>(_data[offset+2]) / 255.f;

          /* RGB -> XYZ easyRGB.com */
          if ( var_R > 0.04045f )
            var_R = pow( ((var_R + 0.055f)/1.055f), 2.4f);
          else
            var_R = var_R / 12.92f;
          if ( var_G > 0.04045f )
            var_G = pow( ((var_G + 0.055f)/1.055f), 2.4f);
          else
            var_G = var_G / 12.92f;
          if ( var_B > 0.04045f )
            var_B = pow( ((var_B + 0.055f)/1.055f), 2.4f);
          else
            var_B = var_B / 12.92f;

          var_R = var_R * 100;
          var_G = var_G * 100;
          var_B = var_B * 100;

          /* Observer. = 2°, Illuminant = D65 */
          float var_X = var_R * 0.4124f + var_G * 0.3576f + var_B * 0.1805f;
          float var_Y = var_R * 0.2126f + var_G * 0.7152f + var_B * 0.0722f;
          float var_Z = var_R * 0.0193f + var_G * 0.1192f + var_B * 0.9505f;

          /* XYZ -> Cie Lab easyRGB.com */
          const float ref_X = 95.047f;
          const float ref_Y = 100.f;
          const float ref_Z = 108.883f;
          var_X = var_X / ref_X;
          var_Y = var_Y / ref_Y;
          var_Z = var_Z / ref_Z;

          if ( var_X > 0.008856f )
            var_X = pow(var_X, 0.333333f);
          else
            var_X = (7.787f * var_X) + (16.f/116.f);
          if ( var_Y > 0.008856f )
            var_Y = pow(var_Y , 0.333333f);
          else
            var_Y = (7.787f * var_Y) + (16.f/116.f);
          if ( var_Z > 0.008856f )
            var_Z = pow(var_Z , 0.333333f);
          else
            var_Z = (7.787f * var_Z) + (16.f/116.f);

          L(i,j) = (116 * var_Y) - 16;
          A(i,j) = 500 * (var_X - var_Y);
          B(i,j) = 200 * (var_Y - var_Z);

          offset += 3;
        }
        offset+=padding;
    }
    file.close();
    free(_data);
  }
  else {
    height = width = 0;
    id = -1;
    l = a = b = NULL;
  }
}

Image::ColorImage::ColorImage(unsigned int height, unsigned int width, int id)
  : height(height) {

  set_width( width);
  set_id( id);
  posix_memalign( (void**)&l, MEMALLIGN, height*_width_*sizeof(float));
  posix_memalign( (void**)&a, MEMALLIGN, height*_width_*sizeof(float));
  posix_memalign( (void**)&b, MEMALLIGN, height*_width_*sizeof(float));
}

Image::ColorImage::ColorImage(Image::ColorImage&& other)
  : height( other.height)
  , width( other.width)
  , _width_( other._width_)
  , id( other.id)
  , l( other.l)
  , a( other.a)
  , b( other.b) {

  other.l = other.a = other.b = NULL;
  other.height = other.width = other._width_ = 0;
  other.id = -1;
}

Image::ColorImage& Image::ColorImage::operator=(Image::ColorImage&& other) {

  if( this != &other) {

    if( l != NULL)
      free( l);
    if( a != NULL)
      free( a);
    if( b != NULL)
      free( b);

    l = other.l;
    a = other.a;
    b = other.b;
    height = other.height;
    width = other.width;
    _width_ = other._width_;
    id = other.id;

    other.l = other.a = other.b = NULL;
    other.height = other.width = other._width_ = 0;
    other.id = -1;
  }

  return *this;
}

Image::ColorImage::ColorImage(const Image::ColorImage& other)
  : height( other.height)
  , width( other.width)
  , _width_( other._width_)
  , id( other.id) {

    unsigned int sz = height * _width_;

  if( other.l != NULL) {
    posix_memalign( (void**)&l, MEMALLIGN, height*_width_*sizeof(float));
    std::copy( other.l, other.l+sz, l);
  }
  else
    l = NULL;

  if( other.a != NULL) {
    posix_memalign( (void**)&a, MEMALLIGN, height*_width_*sizeof(float));
    std::copy( other.a, other.a+sz, a);
  }
  else
    a = NULL;

  if( other.b != NULL) {
    posix_memalign( (void**)&b, MEMALLIGN, height*_width_*sizeof(float));
    std::copy( other.b, other.b+sz, b);
  }
  else
    b = NULL;
}

Image::ColorImage& Image::ColorImage::operator=(const Image::ColorImage& other) {

  if( this != &other) {

    unsigned int sz = other.height * other._width_;
    if( !l || !a || !b || other.width != width || other.height != height) {

      if( l != NULL) {
        free( l);
      }
      if( a != NULL) {
        free( a);
      }
      if( b != NULL) {
        free( b);
      }

      height = other.height;
      width = other.width;
      _width_ = other._width_;
      posix_memalign( (void**)&l, MEMALLIGN, sz*sizeof(float));
      posix_memalign( (void**)&a, MEMALLIGN, sz*sizeof(float));
      posix_memalign( (void**)&b, MEMALLIGN, sz*sizeof(float));
    }
    std::copy( other.l, other.l+sz, l);
    std::copy( other.a, other.a+sz, a);
    std::copy( other.b, other.b+sz, b);
  }

  return *this;
}

Image::ColorImage Image::ColorImage::scale_image(float scale_factor) const {

#if _DEBUG == 1
  assert( scale_factor > 0.f);
#endif
  Image::ColorImage im( static_cast<unsigned int>( floor(height * scale_factor)),
                       static_cast<unsigned int>( floor(width * scale_factor)), id);

  const fp step = 1.0/static_cast<fp>( scale_factor);

  for( unsigned int i = 0; i < im.height; i++) {
    fp y = i*step;
#if INTERPOLATE_FLAG == 1
    unsigned int y0 = static_cast<unsigned int>( floor(y));
    unsigned int y1 = min(y0 + 1, height - 1);
#else
    unsigned int y0 = min( static_cast<unsigned int>( round(y)), height-1);
#endif
//#pragma simd
    for( unsigned int j = 0; j < im.width; j++) {
      fp x = j * step;
#if INTERPOLATE_FLAG == 1
      unsigned int x0 = static_cast<unsigned int>( floor(x));
      unsigned int x1 = min(x0 + 1, width-1);
      fp f00 = L(y0, x0);
      fp f01 = L(y0, x1);
      fp f10 = L(y1, x0);
      fp f11 = L(y1, x1);
      im.L(i,j) = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                      (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
      f00 = A(y0, x0);
      f01 = A(y0, x1);
      f10 = A(y1, x0);
      f11 = A(y1, x1);
      im.A(i,j) = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                      (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
      f00 = B(y0, x0);
      f01 = B(y0, x1);
      f10 = B(y1, x0);
      f11 = B(y1, x1);
      im.B(i,j) = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                      (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
#else
      unsigned int x0 = min( static_cast<unsigned int>( round(x)), width-1);
      im.L(i,j) = L( y0, x0);
      im.A(i,j) = A( y0, x0);
      im.B(i,j) = B( y0, x0);
#endif
    }
  }

  return im;
}

/* rotate relative to image center */
Image::ColorImage Image::ColorImage::rotate_image(float angle) const {

  Image::ColorImage im( height, width, id);

  const fp xc = static_cast<fp>( width/2);
  const fp yc = static_cast<fp>( height/2);

  const fp cos_alpha = cos( Utils::D2R * angle);
  const fp sin_alpha = sin( Utils::D2R * angle);

  for( unsigned int i=0; i < height; ++i) {

    fp yp = static_cast<fp>( i) - yc;
#pragma simd
    for( unsigned int j=0; j < width; ++j) {

      fp xp = static_cast<fp>( j) - xc;

      fp x = ( xp * cos_alpha + yp * sin_alpha) + xc;
      fp y = ( -xp * sin_alpha + yp * cos_alpha) + yc;

      if( (x >= 0) && (x <= (width-1)) && (y >= 0) && (y <= (height-1)) ) {
#if INTERPOLATE_FLAG == 1
        unsigned int y0 = static_cast< unsigned int>( floor(y));
        unsigned int x0 = static_cast< unsigned int>( floor(x));
        unsigned int y1 = y0 + 1;
        unsigned int x1 = x0 + 1;
        fp f00 = static_cast<fp>( L( y0, x0));
        fp f01 = static_cast<fp>( L( y0, x1));
        fp f10 = static_cast<fp>( L( y1, x0));
        fp f11 = static_cast<fp>( L( y1, x1));
        im.L( i, j) = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                          (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
        f00 = static_cast<fp>( A( y0, x0));
        f01 = static_cast<fp>( A( y0, x1));
        f10 = static_cast<fp>( A( y1, x0));
        f11 = static_cast<fp>( A( y1, x1));
        im.A( i, j) = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                          (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
        f00 = static_cast<fp>( B( y0, x0));
        f01 = static_cast<fp>( B( y0, x1));
        f10 = static_cast<fp>( B( y1, x0));
        f11 = static_cast<fp>( B( y1, x1));
        im.B( i, j) = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                          (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
#else
        unsigned int y0 = static_cast< unsigned int>( round(y));
        unsigned int x0 = static_cast< unsigned int>( round(x));
        im.L( i, j) = L( y0, x0);
        im.A( i, j) = A( y0, x0);
        im.B( i, j) = B( y0, x0);
#endif
      }
      else
        im.L( i, j) = im.A( i, j) = im.B( i, j) = 0.f;
    }
  }

  return im;
}

Image::ColorImage::~ColorImage() {

  if( l != NULL) {
    free( l);
  }
  if( a != NULL) {
    free( a);
  }
  if( b != NULL) {
    free( b);
  }
}

void Image::ColorImage::write_image_to_bitmap(const Image::ColorImage& im, const std::string& file_name) {

  HeaderStr header;
  const unsigned int bytes_per_pixel = 3;
  const unsigned int width = im.width;
  const unsigned int height = im.height;
  const unsigned int padding = ((width*bytes_per_pixel) % 4 == 0) ? 0 : 4 - ((width*bytes_per_pixel) % 4);
  unsigned int padded_width = width*bytes_per_pixel + padding;

  header.magic_number = 19778;
  header.size = 14 + 40 + padded_width * height;
  header.reserved = 0;
  header.offset = 54;

  header.dibSize = 40;
  header.width = width;
  header.height = height;
  header.plane = 1;
  header.bit_per_pixel = 24;
  header.compression = 0;
  header.data_size = padded_width * height;
  header.hor_res = 2835; header.vert_res = 2835;
  header.color_number = 0;
  header.important = 0;

  std::ofstream file( file_name.c_str(), ios::out|ios::binary|ios::trunc);

    if (file.is_open())	 {

      file.write( (char*)&header, sizeof(HeaderStr));

      /* writting padded data */
      for(int i=static_cast<int>(height-1);i>=0; i--) {
        for(unsigned int j=0; j < width; j++) {

          /* Cie Lab -> XYZ EasyRGB.com */
          float var_Y = ( im.L(i,j) + 16.f ) / 116.f;
          float var_X = im.A(i,j) / 500.f + var_Y;
          float var_Z = var_Y - im.B(i,j) / 200.f;

          if ( pow(var_Y,3.f) > 0.008856f )
            var_Y = pow(var_Y, 3.f);
          else
            var_Y = ( var_Y - 16.f / 116.f ) / 7.787f;
          if ( pow(var_X,3.f) > 0.008856f )
            var_X = pow(var_X,3.f);
          else
            var_X = ( var_X - 16.f / 116.f ) / 7.787f;
          if ( pow( var_Z, 3.f) > 0.008856f )
            var_Z = pow(var_Z,3.f);
          else
            var_Z = ( var_Z - 16.f / 116.f ) / 7.787f;

          /*Observer= 2°, Illuminant= D65*/
          const float ref_X = 95.047f;
          const float ref_Y = 100.f;
          const float ref_Z = 108.883f;
          var_X = ref_X * var_X;
          var_Y = ref_Y * var_Y;
          var_Z = ref_Z * var_Z;

          /* XYZ -> RGB EasyRGB.com */

          var_X = var_X / 100.f; //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
          var_Y = var_Y / 100.f; //Y from 0 to 100.000
          var_Z = var_Z / 100.f; //Z from 0 to 108.883

          fp var_R = var_X *  3.2406f + var_Y * -1.5372f + var_Z * -0.4986f;
          fp var_G = var_X * -0.9689f + var_Y *  1.8758f + var_Z *  0.0415f;
          fp var_B = var_X *  0.0557f + var_Y * -0.2040f + var_Z *  1.0570f;

          if (var_R > 0.0031308f)
            var_R = 1.055f * pow(var_R, (1.f/2.4f)) - 0.055f;
          else
            var_R = 12.92f * var_R;
          if (var_G > 0.0031308f)
            var_G = 1.055f * pow(var_G, (1.f/2.4f)) - 0.055f;
          else
            var_G = 12.92f * var_G;
          if (var_B > 0.0031308f)
            var_B = 1.055f * pow(var_B, (1.f/2.4f)) - 0.055f;
          else
            var_B = 12.92f * var_B;

          var_R = var_R * 255.f;
          var_R = var_R > 255.f ? 255.f : var_R;
          var_R = var_R < 0.f ? 0.f : var_R;
          var_G = var_G * 255.f;
          var_G = var_G > 255.f ? 255.f : var_G;
          var_G = var_G < 0.f ? 0.f : var_G;
          var_B = var_B * 255.f;
          var_B = var_B > 255.f ? 255.f : var_B;
          var_B = var_B < 0.f ? 0.f : var_B;
          file<< static_cast<unsigned char>( round(var_B));
          file<< static_cast<unsigned char>( round(var_G));
          file<< static_cast<unsigned char>( round(var_R));
        }

        /* padding */
        for(unsigned int k=0;k<padding;k++) {
          file<<0;
        }
      }

      file.close();
    }
}

Image::ColorImage Image::ColorImage::gaussian_smoother(const Image::ColorImage& im) {

  Image::ColorImage smooth(im.height, im.width, im.id);

  static const int filter[5][5] = {
    1,  4,  7,  4, 1,
    4, 16, 26, 16, 4,
    7, 26, 41, 26, 7,
    4, 16, 26, 16, 4,
    1,  4,  7,  4, 1
  };

  static const float filter_sum = 273.f;
  static const int half_filter_width = 2;
  static const int half_filter_height = 2;
  static const float sigma = 1.0f;
  const unsigned int _i = im.height-2;
  const unsigned int _j = im.width-2;

  for( unsigned int i=2; i<_i; ++i) {
#pragma simd
    for( unsigned int j=2; j<_j; ++j) {
      float sum = 0.f;

      sum += im.L(i-2, j-2) * filter[0][0];
      sum += im.L(i-2, j-1) * filter[0][1];
      sum += im.L(i-2, j) * filter[0][2];
      sum += im.L(i-2, j+1) * filter[0][3];
      sum += im.L(i-2, j+2) * filter[0][4];

      sum += im.L(i-1, j-2) * filter[1][0];
      sum += im.L(i-1, j-1) * filter[1][1];
      sum += im.L(i-1, j) * filter[1][2];
      sum += im.L(i-1, j+1) * filter[1][3];
      sum += im.L(i-1, j+2) * filter[1][4];

      sum += im.L(i, j-2) * filter[2][0];
      sum += im.L(i, j-1) * filter[2][1];
      sum += im.L(i, j) * filter[2][2];
      sum += im.L(i, j+1) * filter[2][3];
      sum += im.L(i, j+2) * filter[2][4];

      sum += im.L(i+1, j-2) * filter[3][0];
      sum += im.L(i+1, j-1) * filter[3][1];
      sum += im.L(i+1, j) * filter[3][2];
      sum += im.L(i+1, j+1) * filter[3][3];
      sum += im.L(i+1, j+2) * filter[3][4];

      sum += im.L(i+2, j-2) * filter[4][0];
      sum += im.L(i+2, j-1) * filter[4][1];
      sum += im.L(i+2, j) * filter[4][2];
      sum += im.L(i+2, j+1) * filter[4][3];
      sum += im.L(i+2, j+2) * filter[4][4];

      smooth.L(i, j) = sum / filter_sum;
    }
  }

  for( unsigned int i=0; i<2; ++i) {
#pragma simd
    for( unsigned int j=0; j<(_j+2); ++j) {
      smooth.L(i, j) = im.L(i, j);
    }
  }

  for( unsigned int i=(im.height-2); i<im.height; ++i) {
#pragma simd
    for( unsigned int j=0; j<(_j+2); ++j) {
      smooth.L(i, j) = im.L(i, j);
    }
  }

  for( unsigned int i=2; i<(im.height-2); ++i) {
    for( unsigned int j=0; j<2; ++j) {
      smooth.L(i, j) = im.L(i, j);
    }
    for( unsigned int j=(im.width-2); j<(im.width); ++j) {
      smooth.L(i, j) = im.L(i, j);
    }
  }

  for( unsigned int i=2; i<_i; ++i) {
#pragma simd assert
    for( unsigned int j=2; j<_j; ++j) {
      float sum = 0.f;

      sum += im.A(i-2, j-2) * filter[0][0];
      sum += im.A(i-2, j-1) * filter[0][1];
      sum += im.A(i-2, j) * filter[0][2];
      sum += im.A(i-2, j+1) * filter[0][3];
      sum += im.A(i-2, j+2) * filter[0][4];

      sum += im.A(i-1, j-2) * filter[1][0];
      sum += im.A(i-1, j-1) * filter[1][1];
      sum += im.A(i-1, j) * filter[1][2];
      sum += im.A(i-1, j+1) * filter[1][3];
      sum += im.A(i-1, j+2) * filter[1][4];

      sum += im.A(i, j-2) * filter[2][0];
      sum += im.A(i, j-1) * filter[2][1];
      sum += im.A(i, j) * filter[2][2];
      sum += im.A(i, j+1) * filter[2][3];
      sum += im.A(i, j+2) * filter[2][4];

      sum += im.A(i+1, j-2) * filter[3][0];
      sum += im.A(i+1, j-1) * filter[3][1];
      sum += im.A(i+1, j) * filter[3][2];
      sum += im.A(i+1, j+1) * filter[3][3];
      sum += im.A(i+1, j+2) * filter[3][4];

      sum += im.A(i+2, j-2) * filter[4][0];
      sum += im.A(i+2, j-1) * filter[4][1];
      sum += im.A(i+2, j) * filter[4][2];
      sum += im.A(i+2, j+1) * filter[4][3];
      sum += im.A(i+2, j+2) * filter[4][4];

      smooth.A(i, j) = sum / filter_sum;
    }
  }

  for( unsigned int i=0; i<2; ++i) {
#pragma simd assert
    for( unsigned int j=0; j<(_j+2); ++j) {
      smooth.A(i, j) = im.A(i, j);
    }
  }

  for( unsigned int i=(im.height-2); i<im.height; ++i) {
#pragma simd assert
    for( unsigned int j=0; j<(_j+2); ++j) {
      smooth.A(i, j) = im.A(i, j);
    }
  }

  for( unsigned int i=2; i<(im.height-2); ++i) {
    for( unsigned int j=0; j<2; ++j) {
      smooth.A(i, j) = im.A(i, j);
    }
    for( unsigned int j=(im.width-2); j<(im.width); ++j) {
      smooth.A(i, j) = im.A(i, j);
    }
  }

  for( unsigned int i=2; i<_i; ++i) {
#pragma simd assert
    for( unsigned int j=2; j<_j; ++j) {
      float sum = 0.f;

      sum += im.B(i-2, j-2) * filter[0][0];
      sum += im.B(i-2, j-1) * filter[0][1];
      sum += im.B(i-2, j) * filter[0][2];
      sum += im.B(i-2, j+1) * filter[0][3];
      sum += im.B(i-2, j+2) * filter[0][4];

      sum += im.B(i-1, j-2) * filter[1][0];
      sum += im.B(i-1, j-1) * filter[1][1];
      sum += im.B(i-1, j) * filter[1][2];
      sum += im.B(i-1, j+1) * filter[1][3];
      sum += im.B(i-1, j+2) * filter[1][4];

      sum += im.B(i, j-2) * filter[2][0];
      sum += im.B(i, j-1) * filter[2][1];
      sum += im.B(i, j) * filter[2][2];
      sum += im.B(i, j+1) * filter[2][3];
      sum += im.B(i, j+2) * filter[2][4];

      sum += im.B(i+1, j-2) * filter[3][0];
      sum += im.B(i+1, j-1) * filter[3][1];
      sum += im.B(i+1, j) * filter[3][2];
      sum += im.B(i+1, j+1) * filter[3][3];
      sum += im.B(i+1, j+2) * filter[3][4];

      sum += im.B(i+2, j-2) * filter[4][0];
      sum += im.B(i+2, j-1) * filter[4][1];
      sum += im.B(i+2, j) * filter[4][2];
      sum += im.B(i+2, j+1) * filter[4][3];
      sum += im.B(i+2, j+2) * filter[4][4];

      smooth.B(i, j) = sum / filter_sum;
    }
  }

  for( unsigned int i=0; i<2; ++i) {
#pragma simd assert
    for( unsigned int j=0; j<(_j+2); ++j) {
      smooth.B(i, j) = im.B(i, j);
    }
  }

  for( unsigned int i=(im.height-2); i<im.height; ++i) {
#pragma simd assert
    for( unsigned int j=0; j<(_j+2); ++j) {
      smooth.B(i, j) = im.B(i, j);
    }
  }

  for( unsigned int i=2; i<(im.height-2); ++i) {
    for( unsigned int j=0; j<2; ++j) {
      smooth.B(i, j) = im.B(i, j);
    }
    for( unsigned int j=(im.width-2); j<(im.width); ++j) {
      smooth.B(i, j) = im.B(i, j);
    }
  }

  return smooth;
}

fp Image::ColorImage::bc_invariant_correlation( const Image::ColorImage& main /* to be matched*/, const Image::ColorImage& temp /* matcher */,
                             unsigned int ym/*height in main image*/, unsigned int xm/*width in main image*/,
                             float scale, float angle) {

#if _DEBUG == 1
  assert( scale > 0.f);
  assert( (ym < main.height) && (xm < main.width));
#endif
  /* sums for cross-correlation */
  fp f_sum=0.0, f_sum2=0.0, t_sum=0.0, t_sum2=0.0, ft_sum=0.0;
  fp S_c=0.0;
  unsigned int count = 0;

  const unsigned int height = static_cast<unsigned int>( floor( temp.get_height() * scale));
  const unsigned int width = static_cast<unsigned int>( floor( temp.get_width() * scale));

  const fp step = 1.0/static_cast<fp>( scale);

  const fp xc = static_cast<fp>( temp.get_width()/2);
  const fp yc = static_cast<fp>( temp.get_height()/2);
  const fp cos_alpha = cos( Utils::D2R * angle);
  const fp sin_alpha = sin( Utils::D2R * angle);

  fp y = 0.0;
  for( unsigned int i = 0; i < height; ++i) {

    const fp yp = (y - yc) * static_cast<fp>( scale); /* height in template image */
#if INTERPOLATE_CORR_FLAG == 1
    const unsigned int y0 = static_cast< unsigned int>( floor(y));
    const unsigned int y1 = min( y0 + 1, temp.height - 1);
#endif

//#pragma simd assert reduction(+:count) reduction(+:f_sum) reduction(+:f_sum2) reduction(+:t_sum) reduction(+:t_sum2) reduction(+:ft_sum) reduction(+:S_c)
    for( unsigned int j = 0; j < width; ++j) {

      fp x = j*step;
      fp xp = (x - xc) * static_cast<fp>( scale); /* width in template image */
      /* height in main image */
      const fp ypp = ( xp * sin_alpha + yp * cos_alpha) + static_cast<fp>( ym);
      /* width in main image */
      const fp xpp = ( xp * cos_alpha - yp * sin_alpha) + static_cast<fp>( xm);

      if( (ypp >= 0) && (ypp < (main.get_height()-1)) && (xpp >= 0) && (xpp < (main.get_width()-1))) {

        /* get value of pixel from template image */
#if INTERPOLATE_CORR_FLAG == 1
        const unsigned int x0 = static_cast< unsigned int>( floor(x));
        const unsigned int x1 = min( x0 + 1, temp.width -1);
        fp f00 = static_cast<fp>( temp.L( y0, x0));
        fp f01 = static_cast<fp>( temp.L( y0, x1));
        fp f10 = static_cast<fp>( temp.L( y1, x0));
        fp f11 = static_cast<fp>( temp.L( y1, x1));
        float t_val = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                          (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
        f00 = static_cast<fp>( temp.A( y0, x0));
        f01 = static_cast<fp>( temp.A( y0, x1));
        f10 = static_cast<fp>( temp.A( y1, x0));
        f11 = static_cast<fp>( temp.A( y1, x1));
        float t_val_a = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                            (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
        f00 = static_cast<fp>( temp.B( y0, x0));
        f01 = static_cast<fp>( temp.B( y0, x1));
        f10 = static_cast<fp>( temp.B( y1, x0));
        f11 = static_cast<fp>( temp.B( y1, x1));
        float t_val_b = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                            (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
#else
        float t_val = static_cast<float>( temp.L( static_cast<unsigned int>( min( static_cast<unsigned int>(round(y)), temp.get_height()-1)),
                                                  static_cast<unsigned int>( min( static_cast<unsigned int>(round(x)), temp.get_width()-1))) );
        float t_val_a = static_cast<float>( temp.A( static_cast<unsigned int>( min( static_cast<unsigned int>(round(y)), temp.get_height()-1)),
                                                    static_cast<unsigned int>( min( static_cast<unsigned int>(round(x)), temp.get_width()-1))) );
        float t_val_b = static_cast<float>( temp.B( static_cast<unsigned int>( min( static_cast<unsigned int>(round(y)), temp.get_height()-1)),
                                                    static_cast<unsigned int>( min( static_cast<unsigned int>(round(x)), temp.get_width()-1))) );
#endif

        /* get value of pixel from main image */
#if INTERPOLATE_CORR_FLAG == 1
        const unsigned int _y0 = static_cast< unsigned int>( floor(ypp));
        const unsigned int _y1 = _y0 + 1;
        const unsigned int _x0 = static_cast< unsigned int>( floor(xpp));
        const unsigned int _x1 = _x0 + 1;
        fp _f00 = static_cast<fp>( main.L( _y0, _x0));
        fp _f01 = static_cast<fp>( main.L( _y0, _x1));
        fp _f10 = static_cast<fp>( main.L( _y1, _x0));
        fp _f11 = static_cast<fp>( main.L( _y1, _x1));
        float f_val = static_cast<float>( (ypp-_y0) * ((xpp-_x0)*_f11 + (_x1-xpp)*_f10) +
                                          (_y1-ypp) * ((xpp-_x0)*_f01 + (_x1-xpp)*_f00) );
        _f00 = static_cast<fp>( main.A( _y0, _x0));
        _f01 = static_cast<fp>( main.A( _y0, _x1));
        _f10 = static_cast<fp>( main.A( _y1, _x0));
        _f11 = static_cast<fp>( main.A( _y1, _x1));
        float f_val_a = static_cast<float>( (ypp-_y0) * ((xpp-_x0)*_f11 + (_x1-xpp)*_f10) +
                                            (_y1-ypp) * ((xpp-_x0)*_f01 + (_x1-xpp)*_f00) );
        _f00 = static_cast<fp>( main.B( _y0, _x0));
        _f01 = static_cast<fp>( main.B( _y0, _x1));
        _f10 = static_cast<fp>( main.B( _y1, _x0));
        _f11 = static_cast<fp>( main.B( _y1, _x1));
        float f_val_b = static_cast<float>( (ypp-_y0) * ((xpp-_x0)*_f11 + (_x1-xpp)*_f10) +
                                            (_y1-ypp) * ((xpp-_x0)*_f01 + (_x1-xpp)*_f00) );
#else
        float f_val = static_cast<float>( main.L( static_cast<unsigned int>( static_cast<unsigned int>(round(ypp))),
                                                  static_cast<unsigned int>( static_cast<unsigned int>(round(xpp))) ));
        float f_val_a = static_cast<float>( main.A( static_cast<unsigned int>( static_cast<unsigned int>(round(ypp))),
                                                    static_cast<unsigned int>( static_cast<unsigned int>(round(xpp))) ));
        float f_val_b = static_cast<float>( main.B( static_cast<unsigned int>( static_cast<unsigned int>(round(ypp))),
                                                    static_cast<unsigned int>( static_cast<unsigned int>(round(xpp))) ));
#endif

        /* sums for cross-correlation */
        t_sum += t_val;
        t_sum2 += t_val * t_val;
        f_sum += f_val;
        f_sum2 += f_val * f_val;
        ft_sum += f_val * t_val;

        S_c += sqrt( pow( t_val_a - f_val_a, 2) + pow( t_val_b - f_val_b, 2));

        count++;
      }
    }

    y += step;
  }
  S_c = 1.f - (S_c/(200.f*sqrt(2.f)*count));
  fp S_l = (ft_sum - (f_sum/count*t_sum)) / sqrt( (f_sum2 - (f_sum/count*f_sum)) * (t_sum2 - (t_sum/count*t_sum)));

  if( S_l > 4.f || S_l < 0.f)
    return 0;

  return pow(S_l, _alpha_) * pow(S_c, _beta_);
}


void Image::ColorImage::tag( Image::ColorImage& main, const Image::ColorImage& temp,
                      unsigned int ym, unsigned int xm, float scale, float angle) {

#if _DEBUG == 1
  assert( scale > 0.f);
  assert( (ym < main.height) && (xm < main.width));
#endif

  const unsigned int height = static_cast<unsigned int>( floor( temp.height * scale));
  const unsigned int width = static_cast<unsigned int>( floor( temp.width * scale));

  const fp step = 1.0/static_cast<fp>( scale);

  const fp xc = static_cast<fp>( temp.width/2);
  const fp yc = static_cast<fp>( temp.height/2);
  const fp cos_alpha = cos( Utils::D2R * angle);
  const fp sin_alpha = sin( Utils::D2R * angle);

  fp y = 0.0;
  for( unsigned int i = 0; i < height; ++i) {

    fp yp = (y - yc) * static_cast<fp>( scale); /* height in template image */
#if INTERPOLATE_CORR_FLAG == 1
    const unsigned int y0 = static_cast< unsigned int>( floor(y));
    const unsigned int y1 = min( y0 + 1, temp.height - 1);
#endif

//#pragma simd assert
    for( unsigned int j = 0; j < width; ++j) {

      fp x = j*step;
      fp xp = (x - xc) * static_cast<fp>( scale); /* width in template image */
      /* height in main image */
      const fp ypp = ( xp * sin_alpha + yp * cos_alpha) + static_cast<fp>( ym);
      /* width in main image */
      const fp xpp = ( xp * cos_alpha - yp * sin_alpha) + static_cast<fp>( xm);

      if( (ypp >= 0) && (ypp < (main.height-1)) && (xpp >= 0) && (xpp < (main.width-1))) {
        /********************TAG L value****************************************/
        /* get value of pixel from template image */
#if INTERPOLATE_CORR_FLAG == 1
        const unsigned int x0 = static_cast< unsigned int>( floor(x));
        const unsigned int x1 = x0 + 1;
        fp f00 = static_cast<fp>( temp.L( y0, x0));
        fp f01 = static_cast<fp>( temp.L( y0, x1));
        fp f10 = static_cast<fp>( temp.L( y1, x0));
        fp f11 = static_cast<fp>( temp.L( y1, x1));
        float t_val = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                          (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
#else
        float t_val = static_cast<float>( temp.L( static_cast<unsigned int>( min( static_cast<unsigned int>(round(y)), temp.height-1)),
                                                  static_cast<unsigned int>( min( static_cast<unsigned int>(round(x)), temp.width-1))) );
#endif
        main.L( static_cast< unsigned int>( round( ypp)), static_cast< unsigned int>( round( xpp)) ) = t_val;

        /********************TAG a value****************************************/
#if INTERPOLATE_CORR_FLAG == 1
        f00 = static_cast<fp>( temp.A( y0, x0));
        f01 = static_cast<fp>( temp.A( y0, x1));
        f10 = static_cast<fp>( temp.A( y1, x0));
        f11 = static_cast<fp>( temp.A( y1, x1));
        t_val = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                    (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
#else
        t_val = static_cast<float>( temp.A( static_cast<unsigned int>( min( static_cast<unsigned int>(round(y)), temp.height-1)),
                                            static_cast<unsigned int>( min( static_cast<unsigned int>(round(x)), temp.width-1))) );
#endif
        main.A( static_cast< unsigned int>( round( ypp)), static_cast< unsigned int>( round( xpp)) ) = t_val;

        /********************TAG b value****************************************/
#if INTERPOLATE_CORR_FLAG == 1
        f00 = static_cast<fp>( temp.B( y0, x0));
        f01 = static_cast<fp>( temp.B( y0, x1));
        f10 = static_cast<fp>( temp.B( y1, x0));
        f11 = static_cast<fp>( temp.B( y1, x1));
        t_val = static_cast<float>( (y-y0) * ((x-x0)*f11 + (x1-x)*f10) +
                                    (y1-y) * ((x-x0)*f01 + (x1-x)*f00) );
#else
        t_val = static_cast<float>( temp.B( static_cast<unsigned int>( min( static_cast<unsigned int>(round(y)), temp.height-1)),
                                            static_cast<unsigned int>( min( static_cast<unsigned int>(round(x)), temp.width-1))) );
#endif
        main.B( static_cast< unsigned int>( round( ypp)), static_cast< unsigned int>( round( xpp)) ) = t_val;
      } // if
    } // for width
    y += step;
  }
}


void Image::circle_pix_mean( unsigned int yc, unsigned int xc, unsigned int dx,
                              unsigned int r,
                              const Image::ColorImage& im, fp* _l, fp* _a, fp* _b) {

  unsigned int count = 0;
  float p;
  unsigned int x, y;
  //fp *aux;

#if _DEBUG == 1
  /* return 0 if no complete circle fits the image */
  assert( xc < r);
  assert( (xc+dx+r) >= im.get_width());
  assert( yc < r);
  assert( (yc+r) >= im.get_width());
#endif

  if( r == 0) {
    count = 1;
//    aux = (im.l + yc*im._width_ + xc);
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _l[i] = aux[i];
    _l[0:dx] = (im.l + yc*im._width_)[xc:dx];
//    aux = (im.a + yc*im._width_ + xc);
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _a[i] = aux[i];
    _a[0:dx] = (im.a + yc*im._width_)[xc:dx];
//    aux = (im.b + yc*im._width_ + xc);
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _b[i] = aux[i];
    _b[0:dx] = (im.b + yc*im._width_)[xc:dx];
  }
  else {
    x = 0;
    y = r;
    /* cross-tip points */
//    aux = im.l + (r+yc)*im._width_ + xc;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _l[i] = aux[i];
//    aux = im.l + (-r+yc)*im._width_ + xc;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _l[i] += aux[i];
//    aux = im.l + yc*im._width_ + xc+r;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _l[i] += aux[i];
//    aux = im.l + yc*im._width_ + xc-r;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _l[i] += aux[i];
//    for(unsigned int i=0;i<dx;i++) {
//      _l[i] = (im.l + (r+yc)*im._width_ + xc)[i];
//      _l[i] += (im.l + (-r+yc)*im._width_ + xc)[i];
//      _l[i] += (im.l + yc*im._width_ + xc + r)[i];
//      _l[i] += (im.l + yc*im._width_ + xc - r)[i];
//    }
    _l[0:dx] = (im.l + (r+yc)*im._width_)[xc:dx] + (im.l + (-r+yc)*im._width_)[xc:dx]
               + (im.l + yc*im._width_)[xc+r:dx] + (im.l + yc*im._width_)[xc-r:dx];

//    aux = im.a + (r+yc)*im._width_ + xc;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _a[i] = aux[i];
//    aux = im.a + (-r+yc)*im._width_ + xc;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _a[i] += aux[i];
//    aux = im.a + yc*im._width_ + xc+r;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _a[i] += aux[i];
//    aux = im.a + yc*im._width_ + xc-r;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _a[i] += aux[i];
//    for(unsigned int i=0;i<dx;i++) {
//      _a[i] = (im.a + (r+yc)*im._width_ + xc)[i];
//      _a[i] += (im.a + (-r+yc)*im._width_ + xc)[i];
//      _a[i] += (im.a + yc*im._width_ + xc + r)[i];
//      _a[i] += (im.a + yc*im._width_ + xc - r)[i];
//    }
    _a[0:dx] = (im.a + (r+yc)*im._width_)[xc:dx] + (im.a + (-r+yc)*im._width_)[xc:dx]
               + (im.a + yc*im._width_)[xc+r:dx] + (im.a + yc*im._width_)[xc-r:dx];
//    aux = im.b + (r+yc)*im._width_ + xc;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _b[i] = aux[i];
//    aux = im.b + (-r+yc)*im._width_ + xc;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _b[i] += aux[i];
//    aux = im.b + yc*im._width_ + xc+r;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _b[i] += aux[i];
//    aux = im.b + yc*im._width_ + xc-r;
//#pragma simd assert
//    for(unsigned int i=0;i<dx;i++)
//      _b[i] += aux[i];
//    for(unsigned int i=0;i<dx;i++) {
//      _b[i] = (im.b + (r+yc)*im._width_ + xc)[i];
//      _b[i] += (im.b + (-r+yc)*im._width_ + xc)[i];
//      _b[i] += (im.b + yc*im._width_ + xc + r)[i];
//      _b[i] += (im.b + yc*im._width_ + xc - r)[i];
//    }
    _b[0:dx] = (im.b + (r+yc)*im._width_)[xc:dx] + (im.b + (-r+yc)*im._width_)[xc:dx]
               + (im.b + yc*im._width_)[xc+r:dx] + (im.b + yc*im._width_)[xc-r:dx];
    count = 4;

    p = 1.25f - static_cast<float>( r);

    while(1) {

      p = p + 2.f*static_cast<float>(++x) + 1.f;
      if( p > 0.f)
        p -= 2.f*static_cast<float>(--y);

      if( x >= y) {
        if( x == y) {
//          aux = (im.l + (yc-y)*im._width_ + xc-x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _l[i] += aux[i];
//          aux = (im.l + (yc-y)*im._width_ + xc+x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _l[i] += aux[i];
//          aux = (im.l + (yc+y)*im._width_ + xc+x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _l[i] += aux[i];
//          aux = (im.l + (yc+y)*im._width_ + xc-x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _l[i] += aux[i];

          //for(unsigned int i=0;i<dx;i++) {
          //  _l[i] += (im.l + (yc-y)*im._width_ + xc-x)[i];
          //  _l[i] += (im.l + (yc-y)*im._width_ + xc+x)[i];
          //  _l[i] += (im.l + (yc+y)*im._width_ + xc+x)[i];
          //  _l[i] += (im.l + (yc+y)*im._width_ + xc-x)[i];
          //}
          _l[0:dx] += (im.l + (yc-y)*im._width_)[xc-x:dx] + (im.l + (yc-y)*im._width_)[xc+x:dx]
                     + (im.l + (yc+y)*im._width_)[xc+x:dx] + (im.l + (yc+y)*im._width_)[xc-x:dx];
//          aux = (im.a + (yc-y)*im._width_ + xc-x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _a[i] += aux[i];
//          aux = (im.a + (yc-y)*im._width_ + xc+x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _a[i] += aux[i];
//          aux = (im.a + (yc+y)*im._width_ + xc+x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _a[i] += aux[i];
//          aux = (im.a + (yc+y)*im._width_ + xc-x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _a[i] += aux[i];

          //for(unsigned int i=0;i<dx;i++) {
          //  _a[i] += (im.a + (yc-y)*im._width_ + xc-x)[i];
          //  _a[i] += (im.a + (yc-y)*im._width_ + xc+x)[i];
          //  _a[i] += (im.a + (yc+y)*im._width_ + xc+x)[i];
          //  _a[i] += (im.a + (yc+y)*im._width_ + xc-x)[i];
          //}
          _a[0:dx] += (im.a + (yc-y)*im._width_)[xc-x:dx] + (im.a + (yc-y)*im._width_)[xc+x:dx]
                     + (im.a + (yc+y)*im._width_)[xc+x:dx] + (im.a + (yc+y)*im._width_)[xc-x:dx];

//          aux = (im.b + (yc-y)*im._width_ + xc-x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _b[i] += aux[i];
//          aux = (im.b + (yc-y)*im._width_ + xc+x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _b[i] += aux[i];
//          aux = (im.b + (yc+y)*im._width_ + xc+x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _b[i] += aux[i];
//          aux = (im.b + (yc+y)*im._width_ + xc-x);
//#pragma simd assert
//          for(unsigned int i=0;i<dx;i++)
//            _b[i] += aux[i];
          //for(unsigned int i=0;i<dx;i++) {
          //  _b[i] += (im.b + (yc-y)*im._width_ + xc-x)[i];
          //  _b[i] += (im.b + (yc-y)*im._width_ + xc+x)[i];
          //  _b[i] += (im.b + (yc+y)*im._width_ + xc+x)[i];
          //  _b[i] += (im.b + (yc+y)*im._width_ + xc-x)[i];
          //}
          _b[0:dx] += (im.b + (yc-y)*im._width_)[xc-x:dx] + (im.b + (yc-y)*im._width_)[xc+x:dx]
                     + (im.b + (yc+y)*im._width_)[xc+x:dx] + (im.b + (yc+y)*im._width_)[xc-x:dx];
          count += 4;
        }
        break;
      }

      /*symmetry points in the other seven octants*/
//      aux = (im.l + (yc+y)*im._width_ + xc+x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _l[i] += aux[i];
//      aux = (im.l + (yc+x)*im._width_ + xc+y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _l[i] += aux[i];
//      aux = (im.l + (yc-x)*im._width_ + xc+y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _l[i] += aux[i];
//      aux = (im.l + (yc-y)*im._width_ + xc+x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _l[i] += aux[i];
//      aux = (im.l + (yc-y)*im._width_ + xc-x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _l[i] += aux[i];
//      aux = (im.l + (yc-x)*im._width_ + xc-y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _l[i] += aux[i];
//      aux = (im.l + (yc+x)*im._width_ + xc-y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _l[i] += aux[i];
//      aux = (im.l + (yc+y)*im._width_ + xc-x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _l[i] += aux[i];
      //for(unsigned int i=0;i<dx;i++) {
      //  _l[i] += (im.l + (yc+y)*im._width_ + xc+x)[i];
      //  _l[i] += (im.l + (yc+x)*im._width_ + xc+y)[i];
      //  _l[i] += (im.l + (yc-x)*im._width_ + xc+y)[i];
      //  _l[i] += (im.l + (yc-y)*im._width_ + xc+x)[i];
      //  _l[i] += (im.l + (yc-y)*im._width_ + xc-x)[i];
      //  _l[i] += (im.l + (yc-x)*im._width_ + xc-y)[i];
      //  _l[i] += (im.l + (yc+x)*im._width_ + xc-y)[i];
      //  _l[i] += (im.l + (yc+y)*im._width_ + xc-x)[i];
      //}
      _l[0:dx] += (im.l + (yc+y)*im._width_)[xc+x:dx] + (im.l + (yc+x)*im._width_)[xc+y:dx]
                  + (im.l + (yc-x)*im._width_)[xc+y:dx] + (im.l + (yc-y)*im._width_)[xc+x:dx]
                  + (im.l + (yc-y)*im._width_)[xc-x:dx] + (im.l + (yc-x)*im._width_)[xc-y:dx]
                  + (im.l + (yc+x)*im._width_)[xc-y:dx] + (im.l + (yc+y)*im._width_)[xc-x:dx];
//      aux = (im.a + (yc+y)*im._width_ + xc+x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _a[i] += aux[i];
//      aux = (im.a + (yc+x)*im._width_ + xc+y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _a[i] += aux[i];
//      aux = (im.a + (yc-x)*im._width_ + xc+y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _a[i] += aux[i];
//      aux = (im.a + (yc-y)*im._width_ + xc+x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _a[i] += aux[i];
//      aux = (im.a + (yc-y)*im._width_ + xc-x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _a[i] += aux[i];
//      aux = (im.a + (yc-x)*im._width_ + xc-y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _a[i] += aux[i];
//      aux = (im.a + (yc+x)*im._width_ + xc-y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _a[i] += aux[i];
//      aux = (im.a + (yc+y)*im._width_ + xc-x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _a[i] += aux[i];
      //for(unsigned int i=0;i<dx;i++) {
      //  _a[i] += (im.a + (yc+y)*im._width_ + xc+x)[i];
      //  _a[i] += (im.a + (yc+x)*im._width_ + xc+y)[i];
      //  _a[i] += (im.a + (yc-x)*im._width_ + xc+y)[i];
      //  _a[i] += (im.a + (yc-y)*im._width_ + xc+x)[i];
      //  _a[i] += (im.a + (yc-y)*im._width_ + xc-x)[i];
      //  _a[i] += (im.a + (yc-x)*im._width_ + xc-y)[i];
      //  _a[i] += (im.a + (yc+x)*im._width_ + xc-y)[i];
      //  _a[i] += (im.a + (yc+y)*im._width_ + xc-x)[i];
      //}
      _a[0:dx] += (im.a + (yc+y)*im._width_)[xc+x:dx] + (im.a + (yc+x)*im._width_)[xc+y:dx]
                  + (im.a + (yc-x)*im._width_)[xc+y:dx] + (im.a + (yc-y)*im._width_)[xc+x:dx]
                  + (im.a + (yc-y)*im._width_)[xc-x:dx] + (im.a + (yc-x)*im._width_)[xc-y:dx]
                  + (im.a + (yc+x)*im._width_)[xc-y:dx] + (im.a + (yc+y)*im._width_)[xc-x:dx];
//      aux = (im.b + (yc+y)*im._width_ + xc+x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _b[i] += aux[i];
//      aux = (im.b + (yc+x)*im._width_ + xc+y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _b[i] += aux[i];
//      aux = (im.b + (yc-x)*im._width_ + xc+y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _b[i] += aux[i];
//      aux = (im.b + (yc-y)*im._width_ + xc+x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _b[i] += aux[i];
//      aux = (im.b + (yc-y)*im._width_ + xc-x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _b[i] += aux[i];
//      aux = (im.b + (yc-x)*im._width_ + xc-y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _b[i] += aux[i];
//      aux = (im.b + (yc+x)*im._width_ + xc-y);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _b[i] += aux[i];
//      aux = (im.b + (yc+y)*im._width_ + xc-x);
//#pragma simd assert
//      for(unsigned int i=0;i<dx;i++)
//        _b[i] += aux[i];
      //for(unsigned int i=0;i<dx;i++) {
      //  _b[i] += (im.b + (yc+y)*im._width_ + xc+x)[i];
      //  _b[i] += (im.b + (yc+x)*im._width_ + xc+y)[i];
      //  _b[i] += (im.b + (yc-x)*im._width_ + xc+y)[i];
      //  _b[i] += (im.b + (yc-y)*im._width_ + xc+x)[i];
      //  _b[i] += (im.b + (yc-y)*im._width_ + xc-x)[i];
      //  _b[i] += (im.b + (yc-x)*im._width_ + xc-y)[i];
      //  _b[i] += (im.b + (yc+x)*im._width_ + xc-y)[i];
      //  _b[i] += (im.b + (yc+y)*im._width_ + xc-x)[i];
      //}
      _b[0:dx] += (im.b + (yc+y)*im._width_)[xc+x:dx] + (im.b + (yc+x)*im._width_)[xc+y:dx]
                  + (im.b + (yc-x)*im._width_)[xc+y:dx] + (im.b + (yc-y)*im._width_)[xc+x:dx]
                  + (im.b + (yc-y)*im._width_)[xc-x:dx] + (im.b + (yc-x)*im._width_)[xc-y:dx]
                  + (im.b + (yc+x)*im._width_)[xc-y:dx] + (im.b + (yc+y)*im._width_)[xc-x:dx];
      count += 8;
    }
  }

//#pragma simd assert
//  for(unsigned int i=0;i<dx;i++) {
//    _l[i] /= count;
//  }
  _l[0:dx] /= count;
//#pragma simd assert
//  for(unsigned int i=0;i<dx;i++) {
//    _a[i] /= count;
//  }
  _a[0:dx] /= count;
//#pragma simd assert
//  for(unsigned int i=0;i<dx;i++) {
//    _b[i] /= count;
//  }
  _b[0:dx] /= count;
}

Sampling::CircularSamplingData Image::circular_sampling( const Image::ColorImage& im,
                                                       uint circle_start, uint circle_step_delta) {

  uint radius = im.get_radius();
  uint count = (radius-circle_start)/circle_step_delta + 1;
  Sampling::CircularSamplingData sdata( circle_start, circle_step_delta, count);

  uint r = circle_start;
  fp l, a, b;

  unsigned int yc = im.get_height()/2;
  unsigned int xc = im.get_width()/2;

  for(uint i=0;i<count;i++) {

    Image::circle_pix_mean( yc, xc, 1, r, im, &l, &a, &b);
    sdata.cis_l[i] = l; sdata.cis_l_S += l; sdata.cis_l_S2 += l*l;
    sdata.cis_a[i] = a;
    sdata.cis_b[i] = b;
    r += circle_step_delta;
  }

  return sdata;
}

void Image::line_pix_mean( unsigned int yc, unsigned int xc, float r, float angle, const Image::ColorImage& im,
                           fp* _l, fp* _a, fp* _b) {
#if _DEBUG == 1
  assert( (angle >= 0.f) && (angle <= 360.f));
  assert( r > 0.f);
  assert( ( yc < im.get_height()) && ( xc < im.get_width()));
#endif

  *_l = 0.0; *_a = 0.0; *_b = 0.0;
  unsigned int count = 0;
  int dx, dy, D;
  int xf, yf, x, y;

  dx = static_cast<int>(round(r * cos( angle * Utils::D2R)));
  dy = static_cast<int>(round(r * sin( angle * Utils::D2R)));

  xf = static_cast<int>(std::min( std::max(static_cast<int>(xc) + dx, 0), static_cast<int>(im.get_width())-1));
  yf = static_cast<int>(std::min( std::max(static_cast<int>(yc) + dy, 0), static_cast<int>(im.get_height())-1));

  if( angle >= 0.f && angle < 45.f) {

    D = 2*dy - dx;
    y = yc;

    for( x = xc; (x <= xf) && (y <= yf); ++x) {

      D += 2*dy;
      if( D > 0) {
        *_l += static_cast<fp>(im.L( ++y, x));
        *_a += static_cast<fp>(im.A( ++y, x));
        *_b += static_cast<fp>(im.B( ++y, x));
        D -= 2*dx;
      }
      else {
        *_l += static_cast<fp>(im.L( y, x));
        *_a += static_cast<fp>(im.A( y, x));
        *_b += static_cast<fp>(im.B( y, x));
      }

      count++;
    }

  }
  else if( angle >= 45.f && angle < 90.f) {

    D = 2*dx + dy;
    x = xc;

    for( y = yc; (y <= yf) && (x <= xf); ++y) {

      D -= 2*dx;
      if( D < 0) {
        *_l += static_cast<fp>(im.L( y, ++x));
        *_a += static_cast<fp>(im.A( y, ++x));
        *_b += static_cast<fp>(im.B( y, ++x));
        D += 2*dy;
      }
      else {
        *_l += static_cast<fp>(im.L( y, x));
        *_a += static_cast<fp>(im.A( y, x));
        *_b += static_cast<fp>(im.B( y, x));
      }

      count++;
    }

  }
  else if( angle >= 90.f && angle < 135.f) {

    D = 2*dx + dy;
    x = xc;

    for( y = yc; (y <= yf) && (x >= xf); ++y) {

      D += 2*dx;
      if( D < 0) {
        *_l += static_cast<fp>(im.L( y, --x));
        *_a += static_cast<fp>(im.A( y, --x));
        *_b += static_cast<fp>(im.B( y, --x));
        D += 2*dy;
      }
      else {
        *_l += static_cast<fp>(im.L( y, x));
        *_a += static_cast<fp>(im.A( y, x));
        *_b += static_cast<fp>(im.B( y, x));
      }

      count++;
    }

  }
  else if( angle >= 135.f && angle < 180.f) {

    D = -2*dy - dx;
    y = yc;

    for( x = xc; (x >= xf) && (y <= yf); --x) {

      D -= 2*dy;
      if( D < 0) {
        *_l += static_cast<fp>(im.L( ++y, x));
        *_a += static_cast<fp>(im.A( ++y, x));
        *_b += static_cast<fp>(im.B( ++y, x));
        D -= 2*dx;
      }
      else {
        *_l += static_cast<fp>(im.L( y, x));
        *_a += static_cast<fp>(im.A( y, x));
        *_b += static_cast<fp>(im.B( y, x));
      }

      count++;
    }
  }
  else if( angle >= 180.f && angle < 225.f) {

    D = 2*dy - dx;
    y = yc;

    for( x = xc; (x >= xf) && (y >= yf); --x) {

      D += 2*dy;
      if( D < 0) {
        *_l += static_cast<fp>(im.L( --y, x));
        *_a += static_cast<fp>(im.A( --y, x));
        *_b += static_cast<fp>(im.B( --y, x));
        D -= 2*dx;
      }
      else {
        *_l += static_cast<fp>(im.L( y, x));
        *_a += static_cast<fp>(im.A( y, x));
        *_b += static_cast<fp>(im.B( y, x));
      }

      count++;
    }
  }
  else if( angle >= 225.f && angle < 270.f) {

    D = 2*dx - dy;
    x = xc;

    for( y = yc; (y >= yf) && (x >= xf); --y) {

      D += 2*dx;
      if( D < 0) {
        *_l += static_cast<fp>(im.L( y, --x));
        *_a += static_cast<fp>(im.A( y, --x));
        *_b += static_cast<fp>(im.B( y, --x));
        D -= 2*dy;
      }
      else {
        *_l += static_cast<fp>(im.L( y, x));
        *_a += static_cast<fp>(im.A( y, x));
        *_b += static_cast<fp>(im.B( y, x));
      }

      count++;
    }
  }
  else if( angle >= 270.f && angle < 315.f) {

    D = 2*dx + dy;
    x = xc;

    for( y = yc; (y >= yf) && (x <= xf); --y) {

      D += 2*dx;
      if( D > 0) {
        *_l += static_cast<fp>(im.L( y, ++x));
        *_a += static_cast<fp>(im.A( y, ++x));
        *_b += static_cast<fp>(im.B( y, ++x));
        D += 2*dy;
      }
      else {
        *_l += static_cast<fp>(im.L( y, x));
        *_a += static_cast<fp>(im.A( y, x));
        *_b += static_cast<fp>(im.B( y, x));
      }

      count++;
    }

  }
  else if( angle >= 315.f && angle <= 360.f) {

    D = 2*dy + dx;
    y = yc;

    for( x = xc; (x <= xf) && (y >= yf); ++x) {

      D += 2*dy;
      if( D < 0) {
        *_l += static_cast<fp>(im.L( --y, x));
        *_a += static_cast<fp>(im.A( --y, x));
        *_b += static_cast<fp>(im.B( --y, x));
        D += 2*dx;
      }
      else {
        *_l += static_cast<fp>(im.L( y, x));
        *_a += static_cast<fp>(im.A( y, x));
        *_b += static_cast<fp>(im.B( y, x));
      }

      count++;
    }

  }

  *_l /= count;
  *_a /= count;
  *_b /= count;
}

void Image::radial_sampling( const Image::ColorImage& im, unsigned int yc, unsigned int xc, unsigned int r,
                             float rotation_start, float rotation_step_delta, unsigned int rotation_step_count,
                              fp* _l, fp* _a, fp* _b) {

  for( unsigned int k=0; k < rotation_step_count; k++) {
    Image::line_pix_mean( yc, xc, r, rotation_start + k*rotation_step_delta, im, &_l[k], &_a[k], &_b[k]);
  }
}
