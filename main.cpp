#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <cmath>

#include "sampling.h"
#include "utils.h"
#include "colorimage.h"
#include "ControlDict.h"

#include "omp.h"

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

#if SHOW_FILTERS == 1
  std::vector< std::tuple< unsigned int/*y*/, unsigned int/*x*/> > first_grade_pixels;
  std::vector< std::tuple< unsigned int/*y*/, unsigned int/*x*/> > second_grade_pixels;
#endif

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

  if( parameters.nb_threads >= 1)
    omp_set_num_threads(parameters.nb_threads);

#if _DEBUG == 1
  try {
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
  /* circular sampling templates */
  j = 0;
  for(const Image::ColorImage& temp : templates) {

    Sampling::CircularSamplingData cs = Image::circular_sampling( temp, circle_start, circle_step_delta);
    cs.id = j;
    cs.scale = scaling_start;
    template_cis.push_back( std::move(cs));
    for( float s = scaling_start + scaling_step_delta; s <= scaling_end + 0.0003f; s+= scaling_step_delta) {
      Image::ColorImage scaled = temp.scale_image(s);
      cs = Image::circular_sampling( scaled, circle_start, circle_step_delta);
      cs.id = j;
      cs.scale = s;
      template_cis.push_back( std::move(cs));
    }
    j++;
  }
  std::sort( template_cis.begin(), template_cis.end());

  /* radial sampling templates */
  const float rotation_step_delta = ( rotation_end - rotation_start) / rotation_step_count;
  /* radial sampling data */
  Utils::Array2d<fp> template_ras_l( parameters.template_names.size(), rotation_step_count);
  fp *template_ras_l_S, *template_ras_l_S2;
  posix_memalign( (void**)&template_ras_l_S, MEMALLIGN, rotation_step_count*sizeof(fp));
  posix_memalign( (void**)&template_ras_l_S2, MEMALLIGN, rotation_step_count*sizeof(fp));
  Utils::Array2d<fp> template_ras_a( parameters.template_names.size(), rotation_step_count);
  Utils::Array2d<fp> template_ras_b( parameters.template_names.size(), rotation_step_count);
  j = 0;
  for(const Image::ColorImage& temp : templates) {
    Image::radial_sampling( temp, temp.get_height()/2, temp.get_width()/2, temp.get_radius(),
                            rotation_start, rotation_step_delta, rotation_step_count,
                            template_ras_l.get_row(j), template_ras_a.get_row(j), template_ras_b.get_row(j));
    template_ras_l_S[j] = template_ras_l.reduce_row(j);
    template_ras_l_S2[j] = template_ras_l.reduce_row2(j);
    j++;
  }

  Image::ColorImage main_image( parameters.main_image_name);
  const unsigned int min_radius = templates[0].get_radius();
  const unsigned int lowi = min_radius;
  const unsigned int highi = main_image.get_height() - min_radius;
  const unsigned int lowj = min_radius;
  const unsigned int highj = main_image.get_width() - min_radius;

#pragma omp parallel default(shared)
  {
    fp* buff_l, *buff_a, *buff_b;
    posix_memalign( (void**)&buff_l, MEMALLIGN, (highj-lowj)*sizeof(fp));
    posix_memalign( (void**)&buff_a, MEMALLIGN, (highj-lowj)*sizeof(fp));
    posix_memalign( (void**)&buff_b, MEMALLIGN, (highj-lowj)*sizeof(fp));
    fp* buff_l_S;
    posix_memalign( (void**)&buff_l_S, MEMALLIGN, (highj-lowj)*sizeof(fp));
    fp* buff_l_S2;
    posix_memalign( (void**)&buff_l_S2, MEMALLIGN, (highj-lowj)*sizeof(fp));

    unsigned int max_radius = std::ceil( templates[parameters.template_names.size()-1].get_radius() * parameters.max_scale);
    uint count = (max_radius-circle_start)/circle_step_delta + 1;
    Utils::Array2d<fp> main_cis_l( highj-lowj, count);
    Utils::Array2d<fp> main_cis_a( highj-lowj, count);
    Utils::Array2d<fp> main_cis_b( highj-lowj, count);
    fp* aux;
    posix_memalign( (void**)&aux, MEMALLIGN, std::max( count, rotation_step_count)*sizeof(fp));
    fp* main_ras_l, *main_ras_a, *main_ras_b;
    posix_memalign( (void**)&main_ras_l, MEMALLIGN, rotation_step_count*sizeof(fp));
    posix_memalign( (void**)&main_ras_a, MEMALLIGN, rotation_step_count*sizeof(fp));
    posix_memalign( (void**)&main_ras_b, MEMALLIGN, rotation_step_count*sizeof(fp));
    unsigned int k, r1;

#pragma omp for \
  private(i, j)
    for(i=lowi; i < highi; i++) {

      std::vector< std::tuple< unsigned int /*width coord*/, unsigned int /*temp_id*/, float /*scale*/> > cis_pix;

      k = 0;
      r1 = circle_start;
      for( ; k < template_cis[0].cis_n; k++) {
        Image::circle_pix_mean( i, lowj, highj-lowj, r1, main_image, buff_l, buff_a, buff_b);
        main_cis_l.scatter(k,buff_l,0);
        main_cis_a.scatter(k,buff_a,0);
        main_cis_b.scatter(k,buff_b,0);
        r1 += circle_step_delta;
      }

      for(j=0;j<(highj-lowj);j++) {
        buff_l_S[j] = __sec_reduce_add( main_cis_l.get_row(j)[0:k]);
        buff_l_S2[j] = __sec_reduce_add( pow( main_cis_l.get_row(j)[0:k], 2));
      }

      for(j=0;j<(highj-lowj);j++) {
        fp S_mt = __sec_reduce_add( main_cis_l.get_row(j)[0:k] * template_cis[0].cis_l[0:k]);
        fp S_l = (S_mt - buff_l_S[j]*template_cis[0].cis_l_S/k)
                / sqrt( (template_cis[0].cis_l_S2 - pow( template_cis[0].cis_l_S, 2)/k)
                        * (buff_l_S2[j] - pow( buff_l_S[j], 2)/k) );
        aux[0:k] = pow( main_cis_a.get_row(j)[0:k] - template_cis[0].cis_a[0:k], 2)
                    + pow( main_cis_b.get_row(j)[0:k] - template_cis[0].cis_b[0:k], 2);
        fp S_c = __sec_reduce_add( sqrt( aux[0:k]));
        S_c = 1.f - (S_c/(200.f*sqrt(2.f)*k));
        fp cis_corr = pow(S_l, _alpha_) * pow(S_c, _beta_);
        if( cis_corr > th1) {
          cis_pix.push_back( std::make_tuple( j+lowj, template_cis[0].id, scaling_start));
#if SHOW_FILTERS == 1
#pragma omp critical (first)
{
           first_grade_pixels.push_back( std::make_tuple( i, j+lowj));
}
#endif
        }
      }

      for( unsigned int temp_id = 1; temp_id < template_cis.size(); temp_id++) {

        unsigned int off = circle_step_delta * (template_cis[temp_id].cis_n - template_cis[0].cis_n);

        if((off+lowi)>i || (i+off)>=highi)
          break;

        for( ; k < template_cis[temp_id].cis_n; k++) {

          Image::circle_pix_mean( i, lowj+off, highj-lowj-2*off, r1, main_image, buff_l, buff_a, buff_b);
          main_cis_l.scatter(k, buff_l, off);
          buff_l_S[off:(highj-lowj-2*off)] += buff_l[0:(highj-lowj-2*off)];
          buff_l_S2[off:(highj-lowj-2*off)] += pow( buff_l[0:(highj-lowj-2*off)], 2);
          main_cis_a.scatter(k, buff_a, off);
          main_cis_b.scatter(k, buff_b, off);
          r1 += circle_step_delta;
        }

        for(j=off;j<(highj-lowj-off);j++) {
          fp S_mt = __sec_reduce_add( main_cis_l.get_row(j)[0:k] * template_cis[temp_id].cis_l[0:k]);
          fp S_l = (S_mt - buff_l_S[j]*template_cis[temp_id].cis_l_S/k)
                  / sqrt( (template_cis[temp_id].cis_l_S2 - pow( template_cis[temp_id].cis_l_S, 2)/k)
                          * (buff_l_S2[j] - pow( buff_l_S[j], 2)/k) );
          aux[0:k] = pow( main_cis_a.get_row(j)[0:k] - template_cis[temp_id].cis_a[0:k], 2)
                      + pow( main_cis_b.get_row(j)[0:k] - template_cis[temp_id].cis_b[0:k], 2);
          fp S_c = __sec_reduce_add( sqrt( aux[0:k]));
          S_c = 1.f - (S_c/(200.f*sqrt(2.f)*k));
          fp cis_corr = pow(S_l, _alpha_) * pow(S_c, _beta_);
          if( cis_corr > th1) {
            cis_pix.push_back( std::make_tuple( j+lowj, template_cis[temp_id].id, template_cis[temp_id].scale));
#if SHOW_FILTERS == 1
#pragma omp critical (first)
{
             first_grade_pixels.push_back( std::make_tuple( i, j+lowj));
}
#endif
          }
        }
      } // template_cis


      // Tefi filter using: cis_id, cis_scale, template_ras_l, template_ras_a, template_ras_b
      for( std::vector< std::tuple< unsigned int /* width coord*/, unsigned int /* temp_id */, float /*scale*/> >::iterator it = cis_pix.begin();
           it != cis_pix.end(); ++it) {

        float radius = templates[std::get<1>(*it)].get_radius() * std::get<2>(*it);

        unsigned int jcoord = std::get<0>(*it);
        Image::radial_sampling( main_image, i, jcoord, static_cast<int>(std::floor(radius)) ,
                                rotation_start, rotation_step_delta, rotation_step_count,
                                main_ras_l, main_ras_a, main_ras_b);
        fp S_f = __sec_reduce_add( main_ras_l[0:rotation_step_count]);
        fp S_f2 = __sec_reduce_add( pow( main_ras_l[0:rotation_step_count], 2));

        fp S_t = template_ras_l_S[std::get<1>(*it)];
        fp S_t2 = template_ras_l_S2[std::get<1>(*it)];

        fp* t_ras_l = template_ras_l.get_row( std::get<1>(*it));
        fp* t_ras_a = template_ras_a.get_row( std::get<1>(*it));
        fp* t_ras_b = template_ras_b.get_row( std::get<1>(*it));

        fp angle = rotation_start;
        fp best_angle = rotation_start;
        fp max_ras_corr = 0;
        for( k=0; k<rotation_step_count; k++) {

          fp S_ft = 0;
          S_ft += __sec_reduce_add( t_ras_l[k:(rotation_step_count-k)] * main_ras_l[0:(rotation_step_count-k)]);
          S_ft += __sec_reduce_add( t_ras_l[0:k] * main_ras_l[(rotation_step_count-k):k]);

          fp S_l  =  (S_ft - S_f*S_t/rotation_step_count) /
              sqrt( (S_f2 - S_f*S_f/rotation_step_count) * (S_t2 - S_t*S_t/rotation_step_count) );

          aux[k:(rotation_step_count-k)] = pow( main_ras_a[0:(rotation_step_count-k)] - t_ras_a[k:(rotation_step_count-k)],2)
                                           + pow( main_ras_b[0:(rotation_step_count-k)] - t_ras_b[k:(rotation_step_count-k)],2);
          aux[0:k] = pow( main_ras_a[(rotation_step_count-k):k] - t_ras_a[0:k], 2)
                     + pow( main_ras_b[(rotation_step_count-k):k] - t_ras_b[0:k], 2);
          fp S_c = __sec_reduce_add( sqrt( aux[0:rotation_step_count]));
          S_c = 1.f - (S_c/(200.f*sqrt(2.f)*rotation_step_count));
          fp ras_corr = pow(S_l, _alpha_) * pow(S_c, _beta_);
          if( ras_corr > max_ras_corr) {
            max_ras_corr = ras_corr;
            best_angle = angle;
          }

          angle += rotation_step_delta;
        }


        if( max_ras_corr > th2) {
#if SHOW_FILTERS == 1
#pragma omp critical (second)
{
             second_grade_pixels.push_back( std::make_tuple( i, std::get<0>(*it)));
}
#endif
        }

      }

    } // i

    free(main_ras_l); free(main_ras_a); free(main_ras_b);
    free(buff_l); free(buff_a); free(buff_b);
    free(buff_l_S);
    free(buff_l_S2);
    free(aux);
  }

  free(template_ras_l_S);
  free(template_ras_l_S2);

#if SHOW_FILTERS == 1
  Image::ColorImage mask_image1( main_image);
  for( std::vector< std::tuple< unsigned int, unsigned int> >::iterator it = first_grade_pixels.begin();
       it != first_grade_pixels.end(); ++it) {
    mask_image1.L(std::get<0>(*it), std::get<1>(*it)) = LMAGENTA;
    mask_image1.A(std::get<0>(*it), std::get<1>(*it)) = AMAGENTA;
    mask_image1.B(std::get<0>(*it), std::get<1>(*it)) = BMAGENTA;
  }

  Image::ColorImage::write_image_to_bitmap( mask_image1, "m1_.bmp");

  Image::ColorImage mask_image2( main_image);
  for( std::vector< std::tuple< unsigned int, unsigned int> >::iterator it = second_grade_pixels.begin();
       it != second_grade_pixels.end(); ++it) {
    mask_image2.L(std::get<0>(*it), std::get<1>(*it)) = LMAGENTA;
    mask_image2.A(std::get<0>(*it), std::get<1>(*it)) = AMAGENTA;
    mask_image2.B(std::get<0>(*it), std::get<1>(*it)) = BMAGENTA;
  }

  Image::ColorImage::write_image_to_bitmap( mask_image2, "m2_.bmp");
#endif

#if _DEBUG == 1
  }
  catch (int e) {
    std::cout<< "An exception occurred. Exception Nr. " << e << '\n';
  }
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
