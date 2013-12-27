#ifndef COLORIMAGE_H
#define COLORIMAGE_H

#include "settings.h"
#include "sampling.h"

#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdint.h>
#include <stdlib.h>
#if _DEBUG == 1
#include <assert.h>
#endif

namespace Image {

  typedef struct __attribute__ ((__packed__)){
      //BMP header
      uint16_t magic_number;
      uint32_t size;
      uint32_t reserved;
      uint32_t offset;
      //DIB header
      uint32_t dibSize;
      uint32_t width;
      uint32_t height;
      uint16_t plane;
      uint16_t bit_per_pixel;
      uint32_t compression;
      uint32_t data_size;
      uint32_t hor_res;
      uint32_t vert_res;
      uint32_t color_number;
      uint32_t important;
  }HeaderStr;

  class ColorImage {

    unsigned int height;
    unsigned int width;
    /* actual width, so that rows start alligned */
    unsigned int _width_;
    int id;

    float* l;
    float* a;
    float* b;
  public:
    ColorImage(const std::string file_name);
    ColorImage(unsigned int height, unsigned int width, int id);
    ColorImage(ColorImage&& other);
    ColorImage& operator=(ColorImage&& other);
    ColorImage(const ColorImage& other);
    ColorImage& operator=(const ColorImage& other);
    ColorImage scale_image(float scale_factor) const;
    ColorImage rotate_image(float angle) const;
    ~ColorImage();

    float& L(unsigned int i, unsigned int j);
    float L(unsigned int i, unsigned int j) const;
    float& A(unsigned int i, unsigned int j);
    float A(unsigned int i, unsigned int j) const;
    float& B(unsigned int i, unsigned int j);
    float B(unsigned int i, unsigned int j) const;

    void set_width(unsigned int width);
    unsigned int get_width() const;
    void set_height(unsigned int height);
    unsigned int get_height() const;
    void set_id(unsigned int id);
    unsigned int get_id() const;
    unsigned int get_radius() const;

    static void write_image_to_bitmap(const ColorImage& im, const std::string& file_name);
    static ColorImage gaussian_smoother(const ColorImage& im);


    /* debug auxiliary function
     * places the temp image (scaled and rotated), over the main image
     */
    static void tag( ColorImage& main /* to be tagged */, const ColorImage& temp /* tag */,
                     unsigned int y /* height in main image */, unsigned int x /* width int main image */,
                     float scale, float angle);

    /* last stage filter
     * brute force filter, uses scale and angle information to compute correlation
     */
    static fp bc_invariant_correlation( const ColorImage& main /* to be matched */,
                                        const ColorImage& temp /* matcher */,
                                        unsigned int y /* height in main image */, unsigned int x /* width in main image */,
                                        float scale, float angle);

    bool operator<(const ColorImage& other) const;
    friend void circle_pix_mean( unsigned int yc, unsigned int xc, unsigned int dx,
                                  unsigned int r,
                                  const Image::ColorImage& im, fp* _l, fp* _a, fp* _b);
  };

  void circle_pix_mean( unsigned int yc, unsigned int xc, unsigned int dx,
                        unsigned int r,
                        const Image::ColorImage& im, fp* _l, fp* _a, fp* _b);

  void line_pix_mean( unsigned int yc, unsigned int xc, float r, float angle, const Image::ColorImage& im,
                      fp *  _l, fp * _a, fp * _b);

  Sampling::CircularSamplingData circular_sampling( const Image::ColorImage& im,
                                        uint circle_start, uint circle_step_delta);

  void radial_sampling( const Image::ColorImage& im, unsigned int yc, unsigned int xc, unsigned int r,
                        float rotation_start, float rotation_step_delta, unsigned int rotation_step_count,
                        fp* restrict _l, fp* restrict _a, fp* restrict _b);

  void frame_target( unsigned int y1 /* upper left corner height coord */,
                     unsigned int x1 /* upper left corner width coord */,
                     unsigned int th /* target height */, unsigned int tw /* target width */,
                     float angle /* target rotation angle */, Image::ColorImage& im );
} //namespace

inline
float& Image::ColorImage::L(unsigned int i, unsigned int j) {
  return l[i*_width_+j];
}

inline
float Image::ColorImage::L(unsigned int i, unsigned int j) const {
  return l[i*_width_+j];
}

inline
float& Image::ColorImage::A(unsigned int i, unsigned int j) {
  return a[i*_width_+j];
}

inline
float Image::ColorImage::A(unsigned int i, unsigned int j) const {
  return a[i*_width_+j];
}

inline
float& Image::ColorImage::B(unsigned int i, unsigned int j) {
  return b[i*_width_+j];
}

inline
float Image::ColorImage::B(unsigned int i, unsigned int j) const {
  return b[i*_width_+j];
}

inline
void Image::ColorImage::set_width(unsigned int width) {
  this->width = width;
  /* alligned width */
  this->_width_ = (((width*sizeof(float) + (MEMALLIGN-1))/MEMALLIGN)*MEMALLIGN)/sizeof(float);
}

inline
unsigned int Image::ColorImage::get_width() const {
  return this->width;
}

inline
void Image::ColorImage::set_height(unsigned int height) {
  this->height = height;
}

inline
unsigned int Image::ColorImage::get_height() const {
  return this->height;
}

inline
void Image::ColorImage::set_id(unsigned int id) {
  this->id = id;
}

inline
unsigned int Image::ColorImage::get_id() const {
  return this->id;
}

inline
unsigned int Image::ColorImage::get_radius() const {
  return ((width < height ? width : height) - 1)/2;
}

inline
bool Image::ColorImage::operator<(const Image::ColorImage& other) const {
  return get_radius() < other.get_radius();
}

#endif // COLORIMAGE_H
