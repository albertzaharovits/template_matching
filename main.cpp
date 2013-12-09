#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <cmath>

#include "sampling.h"
#include "utils.h"
#include "colorimage.h"
#include "ControlDict.h"

using namespace std;

/*!
 *\struct Parameters
 *\brief This structure holds the application parameters
 */
typedef struct{
    int nb_threads;
    string main_image_name;
    list<std::string> template_names;
    float max_scale;
}Parameters;

/*!
 * \brief Read the parameters
 * \param argc The number of parameters (including the program call).
 * \param argv The array of parameters.
 * \return The parameters have been read correctly.
 */
bool read_parameters(int argc, char* argv[], Parameters& parameters){
    if(argc < 4) return false;

    parameters.nb_threads = atoi(argv[1]);
    if(parameters.nb_threads < 0) return false;

    parameters.max_scale = atof(argv[2]);
    if(parameters.max_scale <= 0.f) return false;

    parameters.main_image_name = string(argv[3]);

    for(uint i=4; i<argc; i++){
        parameters.template_names.push_back(string(argv[i]));
    }
    return true;
}

int main(int argc, char* argv[]) {

  unsigned int i, j;
  Parameters parameters;
  if(!read_parameters(argc, argv, parameters)){
      cout<<"Wrong number of parameters or invalid parameters..."<<endl;
      cout<<"The program must be called with the following parameters:"<<endl;
      cout<<"\t- num_threads: The number of threads"<<endl;
      cout<<"\t- max_scale: The maximum scale that can be applied to the templates in the main image"<<endl;
      cout<<"\t- main_image: The main image path"<<endl;
      cout<<"\t- t1 ... tn: The list of the template paths. Each template separated by a space"<<endl;
      cout<<endl<<"For example : ./run 4 3 img.bmp template1.bmp template2.bmp"<<endl;
      return -1;
  }

#if SHOW_FILTERS == 1
  std::vector< std::tuple< unsigned int/*y*/, unsigned int/*x*/> > first_grade_pixels;
#endif

  vector< Image::ColorImage> templates;
  /* iterates over the pattern images */
  for(const string& template_name : parameters.template_names) {

    /* read a specific pattern image */
    Image::ColorImage template_image( template_name);
    templates.push_back(template_image);
  }

  std::sort( templates.begin(), templates.end());

  /* extra sampling parameters computation*/
  const float scaling_start = 1.f;
  const uint scaling_step_count = floor( parameters.max_scale - scaling_start)/scaling_step_delta + 1u;
  const float scaling_end = scaling_start + (scaling_step_count - 1u) * scaling_step_delta;

  /* circular sampling data */
  std::vector< Sampling::CircularSamplingData > template_cis;
  template_cis.reserve( parameters.template_names.size() * scaling_step_count);
  /* sampling templates */
  for(const Image::ColorImage& temp : templates) {

    Sampling::CircularSamplingData cs = Image::circle_sampling( temp, circle_start, circle_step_delta);
    cs.id = temp.get_id();
    cs.scale = scaling_start;
    template_cis.push_back( std::move(cs));
    for( float s = scaling_start + scaling_step_delta; s <= scaling_end; s+= scaling_step_delta) {
      Image::ColorImage scaled = temp.scale_image(s);
      cs = Image::circle_sampling( scaled, circle_start, circle_step_delta);
      cs.id = temp.get_id();
      cs.scale = s;
      template_cis.push_back( std::move(cs));
    }
  }

  std::sort( template_cis.begin(), template_cis.end());

  Image::ColorImage main_image( parameters.main_image_name);
  unsigned int min_radius = templates[0].get_radius();
  unsigned int lowi = min_radius;
  unsigned int highi = main_image.get_height() - min_radius;
  unsigned int lowj = min_radius;
  unsigned int highj = main_image.get_width() - min_radius;

  { // forces destructors free memory
    fp* buff_l, *buff_a, *buff_b;
    posix_memalign( (void**)&buff_l, MEMALLIGN, (highj-lowj+1)*sizeof(fp));
    posix_memalign( (void**)&buff_a, MEMALLIGN, (highj-lowj+1)*sizeof(fp));
    posix_memalign( (void**)&buff_b, MEMALLIGN, (highj-lowj+1)*sizeof(fp));
    fp* buff_l_S, *buff_a_S, *buff_b_S;
    posix_memalign( (void**)&buff_l_S, MEMALLIGN, (highj-lowj+1)*sizeof(fp));
    posix_memalign( (void**)&buff_a_S, MEMALLIGN, (highj-lowj+1)*sizeof(fp));
    posix_memalign( (void**)&buff_b_S, MEMALLIGN, (highj-lowj+1)*sizeof(fp));
    fp* buff_l_S2, *buff_a_S2, *buff_b_S2;
    posix_memalign( (void**)&buff_l_S2, MEMALLIGN, (highj-lowj+1)*sizeof(fp));
    posix_memalign( (void**)&buff_a_S2, MEMALLIGN, (highj-lowj+1)*sizeof(fp));
    posix_memalign( (void**)&buff_b_S2, MEMALLIGN, (highj-lowj+1)*sizeof(fp));
    fp* cis_corr;
    posix_memalign( (void**)&cis_corr, MEMALLIGN, (highj-lowj+1)*sizeof(fp));

    unsigned int max_radius = std::floor( templates[parameters.template_names.size()-1].get_radius() * parameters.max_scale);
    uint count = (max_radius-circle_start)/circle_step_delta + 1;
    Utils::Array2d<fp> main_l( highj-lowj+1, count);
    Utils::Array2d<fp> main_a( highj-lowj+1, count);
    Utils::Array2d<fp> main_b( highj-lowj+1, count);
    fp* aux;
    posix_memalign( (void**)&aux, MEMALLIGN, count*sizeof(fp));

    for(i=lowi; i < highi; i++) {
      unsigned int k;
      unsigned int r1 = circle_start;
      for( k=0; k < template_cis[0].cis_n; k++) {
        Image::circle_pix_mean( i, lowj, highj-lowj, r1, main_image, buff_l, buff_a, buff_b);
        main_l.scatter(k,buff_l,0);
        main_a.scatter(k,buff_a,0);
        main_b.scatter(k,buff_b,0);
        r1 += circle_step_delta;
      }

      for(j=0;j<(highj-lowj);j++) {
        buff_l_S[j] = main_l.reduce_row(j);
        buff_l_S2[j] = main_l.reduce_row2(j);
      }
      for(j=0;j<(highj-lowj);j++) {
        buff_a_S[j] = main_a.reduce_row(j);
        buff_a_S2[j] = main_a.reduce_row2(j);
      }
      for(j=0;j<(highj-lowj);j++) {
        buff_b_S[j] = main_b.reduce_row(j);
        buff_b_S2[j] = main_b.reduce_row2(j);
      }

      for(j=0;j<(highj-lowj);j++) {
        fp S_mt = __sec_reduce_add( main_l.get_row(j)[0:k] * template_cis[0].cis_l[0:k]);
        fp S_l = (S_mt - buff_l_S[j]*template_cis[0].cis_l_S/k)
                / sqrt( (template_cis[0].cis_l_S2 - pow( template_cis[0].cis_l_S, 2)/k)
                        * (buff_l_S2[j] - pow( buff_l_S[j], 2)/k) );
        aux[0:k] = (pow( main_a.get_row(j)[0:k] - template_cis[0].cis_a[0:k], 2)
                    + pow( main_b.get_row(j)[0:k] - template_cis[0].cis_b[0:k], 2));
        fp S_c = __sec_reduce_add( sqrt( aux[0:k]));
        S_c = 1.f - (S_c/(200.f*sqrt(2.f)*k));
        cis_corr[j] = pow(S_l, _alpha_) * pow(S_c, _beta_);
      }

#if SHOW_FILTERS == 1
      for(j=0;j<(highj-lowj);j++) {
        if( cis_corr[j] > th1)
           first_grade_pixels.push_back( std::make_tuple( i, j+lowj));
      }
#endif

//      for( unsigned int temp_id = 1; temp_id < template_cis.size(); temp_id++) {
//        const Sampling::CircularSamplingData& temp_cis = template_cis[temp_id];
//
//        unsigned int off = circle_step_delta * (template_cis[temp_id-1].cis_n - template_cis[temp_id].cis_n);
//        unsigned int r2 = r1 + off;
//
//        for( ; k < template_cis[temp_id].cis_n; k++) {
//          Image::circle_pix_mean( i, lowj-off, highj-lowj-2*off, r1, main_image, buff_l, buff_a, buff_b);
//          main_l.scatter(k, buff_l, off);
//          buff_l_S[off:(highj-lowj-2*off+1)] += buff_l[0:(highj-lowj-2*off+1)];
//          main_a.scatter(k, buff_a, off);
//          main_b.scatter(k, buff_b, off);
//          r1 += circle_step_delta;
//        }
//      }
//
//      unsigned int off = circle_step_delta;
//      for(unsigned int r=(min_radius+circle_step_delta); r <= std::min(i, max_radius); r+=circle_step_delta) {
//        Image::circle_pix_mean( i, lowj+off, highj-lowj-2*off, r, main_image, buff_l, buff_a, buff_b);
//        main_l.scatter(k,buff_l,off);
//        for(j=off;j<(highj-lowj-2*off+1);j++) {
//          buff_l_S[j] += buff_l[j-off];
//          buff_l_S2[j] += buff_l[j-off]*buff_l[j-off];
//        }
//        main_a.scatter(k,buff_a,off);
//        for(j=off;j<(highj-lowj-2*off+1);j++) {
//          buff_a_S[j] += buff_a[j-off];
//          buff_a_S2[j] += buff_a[j-off]*buff_a[j-off];
//        }
//        main_b.scatter(k,buff_b,off);
//        for(j=off;j<(highj-lowj-2*off+1);j++) {
//          buff_b_S[j] += buff_b[j-off];
//          buff_b_S2[j] += buff_b[j-off]*buff_a[j-off];
//        }
//        k++;
//      }

    }

    free(buff_l); free(buff_a); free(buff_b);
    free(buff_l_S); free(buff_a_S); free(buff_b_S);
    free(buff_l_S2); free(buff_a_S2); free(buff_b_S2);
    free(cis_corr);
    free(aux);
  }

#if SHOW_FILTERS == 1
  Image::ColorImage mask_image1( main_image);
  for( std::vector< std::tuple< unsigned int, unsigned int> >::iterator it = first_grade_pixels.begin();
       it != first_grade_pixels.end(); ++it) {
    mask_image1.L(std::get<0>(*it), std::get<1>(*it)) = LMAGENTA;
    mask_image1.A(std::get<0>(*it), std::get<1>(*it)) = AMAGENTA;
    mask_image1.B(std::get<0>(*it), std::get<1>(*it)) = BMAGENTA;
  }

  Image::ColorImage::write_image_to_bitmap( mask_image1, "m1_.bmp");
#endif

  return 0;
}

//    Image::ColorImage::write_image_to_bitmap(img, "/home/albert/workspace/ayc/original.bmp");
//    Image::ColorImage blur = Image::ColorImage::gaussian_smoother( img);
//    Image::ColorImage::write_image_to_bitmap(blur, "/home/albert/workspace/ayc/blur.bmp");
//    Image::ColorImage::write_image_to_bitmap(img.scale_image(1.2f), "/home/albert/workspace/ayc/scaled12.bmp");
//    Image::ColorImage::write_image_to_bitmap(img.scale_image(0.5f), "/home/albert/workspace/ayc/scaled05.bmp");
//    Image::ColorImage::write_image_to_bitmap(img.rotate_image(10.f), "/home/albert/workspace/ayc/rotate10.bmp");
//    Image::ColorImage::write_image_to_bitmap(img.rotate_image(20.f), "/home/albert/workspace/ayc/rotate20.bmp");

//    Image::ColorImage::write_image_to_bitmap(img.scale_image(1.2f), "scaled12.bmp");
//    Image::ColorImage::write_image_to_bitmap(img.scale_image(0.5f), "scaled05.bmp");
//    Image::ColorImage temp("test_case_2/001template.bmp");
//    Image::ColorImage::tag( img, temp, 500, 500, 2.f, 0.f);
//    fp corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 2.f, 0.f);
//    cout << corr << endl;
//
//    cout<<endl;
//
//    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.9f, 0.f);
//    cout << corr << endl;
//    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.8f, 0.f);
//    cout << corr << endl;
//    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.7f, 0.f);
//    cout << corr << endl;
//
//    cout<<endl;
//
//    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 498, 499, 2.f, 0.f);
//    cout << corr << endl;
//
//    cout<<endl;
//
//    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.7f, 15.f);
//    cout << corr << endl;
//
//    Image::ColorImage::write_image_to_bitmap(img, "tagged.bmp");
