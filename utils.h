#ifndef UTILS_H
#define UTILS_H

#include "settings.h"
#include <cstring>

#define LMAGENTA 60.31993f
#define AMAGENTA 98.25421f
#define BMAGENTA -60.84298f

namespace Utils {

  static const fp PI = 3.14159265;
  static const fp D2R = PI/180.0;

  template< typename T>
  class Array2d {

  private:
    unsigned int height;
    unsigned int width;
    // aligned width
    unsigned int _width_;
    T* data;

  public:
    Array2d( unsigned int height=0, unsigned int width=0);
    ~Array2d();
    Array2d( Array2d<T>&& other);
    Array2d<T>& operator=( Array2d<T>&& other);
    Array2d( const Array2d<T>& other);
    Array2d<T>& operator=( const Array2d<T>& other);

    /* no bounds check */
    T& operator()( unsigned int i, unsigned int j);
    T operator()( unsigned int i, unsigned int j) const;

    Array2d<T>& add( const Array2d<T>& other);
    Array2d<T>& add_2( const Array2d<T>& other);

    /* no bounds check */
    void scatter( unsigned int column, const T* arr, unsigned int offset=0);
    void gather( unsigned int column, T* arr, unsigned int offset=0) const;

    T reduce_row( unsigned int row) const;
    T reduce_row2( unsigned int row) const;

    unsigned int get_width() const;
    unsigned int get_height() const;
    void set_width(unsigned int width);
    void set_height(unsigned int height);

    const T* get_data() const;
    const T* get_row(unsigned int row) const;
    T* get_data();
    T* get_row(unsigned int row);
  };

}

template< typename T>
inline
const T* Utils::Array2d<T>::get_row(unsigned int row) const {
  return data + row*_width_;
}

template< typename T>
inline
T* Utils::Array2d<T>::get_row(unsigned int row) {
  return data + row*_width_;
}

template< typename T>
inline
const T* Utils::Array2d<T>::get_data() const {
  return data;
}

template< typename T>
inline
T* Utils::Array2d<T>::get_data() {
  return data;
}

template< typename T>
Utils::Array2d<T>::Array2d( unsigned int height, unsigned int width)
  : height( height)
  , width( width) {

  _width_ = (((width*sizeof(T) + (MEMALLIGN-1))/MEMALLIGN)*MEMALLIGN)/sizeof(T);
  posix_memalign( (void**)&data, MEMALLIGN, height*_width_*sizeof(T));
  std::memset( data, 0, height*_width_*sizeof(T));
}

template< typename T>
Utils::Array2d<T>::~Array2d() {

  if(data)
    free(data);
}

template< typename T>
Utils::Array2d<T>::Array2d( Utils::Array2d<T>&& other)
  :height( other.height)
  ,width( other.width)
  ,_width_( other._width_)
  ,data( other.data) {

  other_height = other.width = other._width_ = 0;
  other.data = NULL;
}

template< typename T>
Utils::Array2d<T>& Utils::Array2d<T>::operator=( Utils::Array2d<T>&& other) {

  if( this != &other) {

    if(data)
      free(data);

    height = other.height;
    width = other.width;
    _width_ = other._width_;
    data = other.data;

    other.height = other.width = other._width_ = 0;
    other.data = NULL;
  }

  return *this;
}

template< typename T>
Utils::Array2d<T>::Array2d( const Utils::Array2d<T>& other)
  : height( other.height)
  , width( other.width)
  , _width_( other._width_)
  , data( NULL) {

  if( other.data != NULL) {
    posix_memalign( (void**)&data, MEMALLIGN, height*_width_*sizeof(T));
    std::copy(other.data, other.data+(height*_width_), data);
  }
}

template< typename T>
Utils::Array2d<T>& Utils::Array2d<T>::operator=( const Utils::Array2d<T>& other) {

  if( this != &other) {

    if( !data || height != other.height || width != other.width) {

      if( data)
        free( data);

      height = other.height;
      width = other.width;
      _width_ = other._width_;
      posix_memalign( (void**)&data, MEMALLIGN, height*_width_*sizeof(T));
    }

    std::copy(other.data, other.data+(height*_width_), data);
  }

  return *this;
}

template< typename T>
Utils::Array2d<T>& Utils::Array2d<T>::add( const Utils::Array2d<T>& other) {

  if( this != &other) {
#ifdef _DEBUG
    assert( (_height == other._height) && (_width == other._width));
#endif

    unsigned int sz = height * _width_;
    data[0:sz] += other.data[0:sz];
  }

  return *this;
}

template< typename T>
Utils::Array2d<T>& Utils::Array2d<T>::add_2( const Utils::Array2d<T>& other) {

  if( this != &other) {
#ifdef _DEBUG
    assert( (_height == other._height) && (_width == other._width));
#endif

    unsigned int sz = height * _width_;
    data[0:sz] += (other.data[0:sz]*other.data[0:sz]);
  }

  return *this;
}

template< typename T>
void Utils::Array2d<T>::scatter( unsigned int column, const T* arr, unsigned int offset) {
  for( unsigned int i=offset; i < (height-offset); i++) {
    data[i*_width_ + column] = arr[i-offset];
  }
}

template< typename T>
void Utils::Array2d<T>::gather( unsigned int column, T* arr, unsigned int offset) const {
  for( unsigned int i=offset; i < (height-offset); i++) {
    arr[i-offset] = data[i*_width_ + column];
  }
}

template< typename T>
inline
T Utils::Array2d<T>::reduce_row( unsigned int row) const {
  return __sec_reduce_add( (data+row*_width_)[0:_width_] );
}

template< typename T>
inline
T Utils::Array2d<T>::reduce_row2( unsigned int row) const {
  return __sec_reduce_add( pow( (data+row*_width_)[0:_width_], 2) );
}

template< typename T>
inline
T& Utils::Array2d<T>::operator()( unsigned int i, unsigned int j) { return data[i*_width_+j]; }

template< typename T>
inline
T Utils::Array2d<T>::operator()( unsigned int i, unsigned int j) const { return data[i*_width_+j]; }

template< typename T>
inline
unsigned int Utils::Array2d<T>::get_height() const { return height; }

template< typename T>
inline
unsigned int Utils::Array2d<T>::get_width() const { return width; }

template< typename T>
inline
void Utils::Array2d<T>::set_height(unsigned int height) { height = height; }

template< typename T>
inline
void Utils::Array2d<T>::set_width(unsigned int width) {
  width = width;
  _width_ = (((width*sizeof(T) + (MEMALLIGN-1))/MEMALLIGN)*MEMALLIGN)/sizeof(T);
}


#endif // UTILS_H
