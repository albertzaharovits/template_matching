#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cmath>

#include "sampling.h"
#include "utils.h"
#include "colorimage.h"
#include "ControlDict.h"

#include <omp.h>

using namespace std;


typedef DisjointSet::DsCell< std::tuple< int /*template id*/, int /*x*/, int /*y*/, fp /*corr*/,
                                         unsigned int /* width */, unsigned int /* height */, float /* angle */> >* pDsCell;

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

#if SHOW_FILTERS == 1
  std::vector< std::tuple< unsigned int/*y*/, unsigned int/*x*/> > first_grade_pixels;
  std::vector< std::tuple< unsigned int/*y*/, unsigned int/*x*/> > second_grade_pixels;
  std::vector< std::tuple< unsigned int/*y*/, unsigned int/*x*/> > third_grade_pixels;
#endif

  std::vector< pDsCell > results;

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

  vector< Image::ColorImage> templates;
  /* iterates over the pattern images */
  for(const string& template_name : parameters.template_names) {

    /* read a specific pattern image */
    Image::ColorImage template_image( template_name);
#if GAUSSIAN_FLAG == 1
    templates.push_back(Image::ColorImage::gaussian_smoother(std::move(template_image)));
#else
    templates.push_back(std::move(template_image));
#endif
  }
  std::stable_sort( templates.begin(), templates.end());

  /* extra sampling parameters computation*/
  const uint scaling_step_count = floor( parameters.max_scale - scaling_start)/scaling_step_delta + 1u;
  const float scaling_end = scaling_start + (scaling_step_count - 1u) * scaling_step_delta;

  /* adjusting radial sampling step delta, accounting for min template radius */
  circle_step_delta = std::max( circle_step_delta,
                                static_cast<uint>(( std::round(templates[0].get_radius() * scaling_start / start_circle_count))));

  /* circular sampling data */
  std::vector< Sampling::CircularSamplingData > template_cis;
  template_cis.reserve( parameters.template_names.size() * scaling_step_count);
  /* circular sampling templates */
  j = 0;
  for(const Image::ColorImage& temp : templates) {
    for( unsigned int k=0; k < scaling_step_count; k++) {
      float s = scaling_start + k*scaling_step_delta;
      Image::ColorImage scaled = temp.scale_image(s);
      Sampling::CircularSamplingData cs = Image::circular_sampling( scaled, circle_start, circle_step_delta);
      cs.id = j;
      cs.scale = s;
      template_cis.push_back( std::move(cs));
    }
    j++;
  }
  std::stable_sort( template_cis.begin(), template_cis.end());

  /* radial sampling templates */
  const float rotation_step_delta = ( rotation_end - rotation_start) / rotation_step_count;

  /* radial sampling data */
  Utils::Array2d<fp> template_ras_l( parameters.template_names.size(), rotation_step_count);
  fp *template_ras_l_S, *template_ras_l_S2;
  posix_memalign( (void**)&template_ras_l_S, MEMALLIGN, parameters.template_names.size()*sizeof(fp));
  posix_memalign( (void**)&template_ras_l_S2, MEMALLIGN, parameters.template_names.size()*sizeof(fp));
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
#if GAUSSIAN_FLAG == 1
  main_image = std::move( Image::ColorImage::gaussian_smoother(main_image));
#endif
  const unsigned int min_radius = template_cis[0].cis_n * circle_step_delta;
  const unsigned int lowi = min_radius;
  const unsigned int highi = main_image.get_height() - min_radius;
  const unsigned int lowj = min_radius;
  const unsigned int highj = main_image.get_width() - min_radius;

#pragma omp parallel default(shared)
  {
    fp* buff_l;
    posix_memalign( (void**)&buff_l, MEMALLIGN, (highj-lowj)*sizeof(fp));
    fp* buff_l_S;
    posix_memalign( (void**)&buff_l_S, MEMALLIGN, (highj-lowj)*sizeof(fp));
    fp* buff_l_S2;
    posix_memalign( (void**)&buff_l_S2, MEMALLIGN, (highj-lowj)*sizeof(fp));

    unsigned int max_radius = std::round( templates[parameters.template_names.size()-1].get_radius() * parameters.max_scale);
    uint count = (max_radius-circle_start)/circle_step_delta + 1 + 1;
    Utils::Array2d<fp> main_cis_l( highj-lowj, count);
    fp* aux, *aux2;
    posix_memalign( (void**)&aux, MEMALLIGN, std::max( count, rotation_step_count)*sizeof(fp));
    posix_memalign( (void**)&aux2, MEMALLIGN, rotation_step_count*sizeof(fp));
    fp* main_ras_l, *main_ras_a, *main_ras_b;
    posix_memalign( (void**)&main_ras_l, MEMALLIGN, rotation_step_count*sizeof(fp));
    posix_memalign( (void**)&main_ras_a, MEMALLIGN, rotation_step_count*sizeof(fp));
    posix_memalign( (void**)&main_ras_b, MEMALLIGN, rotation_step_count*sizeof(fp));
    unsigned int k, r1;

#pragma omp for \
  private(i, j) schedule(guided)
    for(i=lowi; i < highi; i++) {

      std::vector< std::tuple< unsigned int /*width coord*/, unsigned int /*temp_id*/, float /*scale*/> > cis_pix;

      k = 0;
      r1 = circle_start;
      for( ; k < template_cis[0].cis_n; k++) {
        Image::circle_pix_mean( i, lowj, highj-lowj, r1, main_image, buff_l);
        main_cis_l.scatter(k,buff_l,0);
        r1 += circle_step_delta;
      }

      for(j=0;j<(highj-lowj);j++) {
        fp *m_cis_l = main_cis_l.get_row(j);
        buff_l_S[j] = __sec_reduce_add( m_cis_l[0:k]);
        buff_l_S2[j] = __sec_reduce_add( pow( m_cis_l[0:k], 2));
      }

      for(j=0;j<(highj-lowj);j++) {
        fp *m_cis_l = main_cis_l.get_row(j);
        fp *t_cis_l = template_cis[0].cis_l;
        fp S_mt = __sec_reduce_add( m_cis_l[0:k] * t_cis_l[0:k]);
        fp S_l = (S_mt - buff_l_S[j]*template_cis[0].cis_l_S/k)
                / sqrt( (template_cis[0].cis_l_S2 - pow( template_cis[0].cis_l_S, 2)/k)
                        * (buff_l_S2[j] - pow( buff_l_S[j], 2)/k) );
        fp cis_corr;
        if( S_l > 2.f || S_l < 0.f) {
          cis_corr = 0.f;
        }
        else {
          cis_corr = S_l;
        }
        if( cis_corr > th1) {
          cis_pix.push_back( std::make_tuple( j+lowj, template_cis[0].id, template_cis[0].scale));
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

          Image::circle_pix_mean( i, lowj+off, highj-lowj-2*off, r1, main_image, buff_l);
          main_cis_l.scatter(k, buff_l, off);
          buff_l_S[off:(highj-lowj-2*off)] += buff_l[0:(highj-lowj-2*off)];
          buff_l_S2[off:(highj-lowj-2*off)] += pow( buff_l[0:(highj-lowj-2*off)], 2);
          r1 += circle_step_delta;
        }

        for(j=off;j<(highj-lowj-off);j++) {
          fp *m_cis_l = main_cis_l.get_row(j);
          fp *t_cis_l = template_cis[temp_id].cis_l;
          fp S_mt = __sec_reduce_add( m_cis_l[0:k] * t_cis_l[0:k]);
          fp S_l = (S_mt - buff_l_S[j]*template_cis[temp_id].cis_l_S/k)
                  / sqrt( (template_cis[temp_id].cis_l_S2 - pow( template_cis[temp_id].cis_l_S, 2)/k)
                          * (buff_l_S2[j] - pow( buff_l_S[j], 2)/k) );
          fp cis_corr;
          if( S_l > 2.f || S_l < 0.f) {
            cis_corr = 0.f;
          }
          else {
            cis_corr = S_l;
          }
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


      for( std::vector< std::tuple< unsigned int /* width coord*/, unsigned int /* temp_id */, float /*scale*/> >::iterator it = cis_pix.begin();
           it != cis_pix.end(); ++it) {

        float radius = templates[std::get<1>(*it)].get_radius() * std::get<2>(*it);

        Image::radial_sampling( main_image, i, std::get<0>(*it), static_cast<int>(std::floor(radius)) ,
                                rotation_start, rotation_step_delta, rotation_step_count,
                                main_ras_l, main_ras_a, main_ras_b);
        fp S_f = __sec_reduce_add( main_ras_l[0:rotation_step_count]);
        fp S_f2 = __sec_reduce_add( pow( main_ras_l[0:rotation_step_count], 2));

        fp S_t = template_ras_l_S[std::get<1>(*it)];
        fp S_t2 = template_ras_l_S2[std::get<1>(*it)];

        fp* t_ras_l = template_ras_l.get_row( std::get<1>(*it));
        fp* t_ras_a = template_ras_a.get_row( std::get<1>(*it));
        fp* t_ras_b = template_ras_b.get_row( std::get<1>(*it));

        for( k=0; k<rotation_step_count; k++) {

          fp S_ft = 0;
          S_ft += __sec_reduce_add( t_ras_l[0:(rotation_step_count-k)] * main_ras_l[k:(rotation_step_count-k)]);
          S_ft += __sec_reduce_add( t_ras_l[(rotation_step_count-k):k] * main_ras_l[0:k]);

          fp S_l  =  (S_ft - S_f*S_t/rotation_step_count) /
              sqrt( (S_f2 - S_f*S_f/rotation_step_count) * (S_t2 - S_t*S_t/rotation_step_count) );

          aux[k:(rotation_step_count-k)] = pow( main_ras_a[k:(rotation_step_count-k)] - t_ras_a[0:(rotation_step_count-k)],2)
                                           + pow( main_ras_b[k:(rotation_step_count-k)] - t_ras_b[0:(rotation_step_count-k)],2);
          aux[0:k] = pow( main_ras_a[0:k] - t_ras_a[(rotation_step_count-k):k], 2)
                     + pow( main_ras_b[0:k] - t_ras_b[(rotation_step_count-k):k], 2);
          fp S_c = __sec_reduce_add( sqrt( aux[0:rotation_step_count]));
          S_c = 1.f - (S_c/(200.f*sqrt(2.f)*rotation_step_count));
          aux2[k] = pow(S_l, _alpha_) * pow(S_c, _beta_);
        }
        unsigned int maxi = __sec_reduce_max_ind( aux2[0:rotation_step_count]);

        if( aux2[maxi] < th2)
          continue;

        fp angle = rotation_start + maxi * rotation_step_delta;

#if SHOW_FILTERS == 1
#pragma omp critical (second)
{
        second_grade_pixels.push_back( std::make_tuple( i, std::get<0>(*it)));
}
#endif
        fp best_scale = std::get<2>(*it);
        fp best_angle = angle;

        fp corr = Image::ColorImage::bc_invariant_correlation( main_image, templates[std::get<1>(*it)],
                                       i, std::get<0>(*it), std::get<2>(*it), angle);

        fp corr2 = Image::ColorImage::bc_invariant_correlation( main_image, templates[std::get<1>(*it)],
                                       i, std::get<0>(*it), std::get<2>(*it), angle+(rotation_step_delta/2));
        if( corr2 > corr) {
          best_angle = angle+(rotation_step_delta/2);
          best_scale = std::get<2>(*it);
          corr = corr2;
        }

        corr2 = Image::ColorImage::bc_invariant_correlation( main_image, templates[std::get<1>(*it)],
                                     i, std::get<0>(*it), std::get<2>(*it), angle-(rotation_step_delta/2));
        if( corr2 > corr) {
          best_angle = angle-(rotation_step_delta/2);
          best_scale = std::get<2>(*it);
          corr = corr2;
        }

        corr2 = Image::ColorImage::bc_invariant_correlation( main_image, templates[std::get<1>(*it)],
                                       i, std::get<0>(*it), std::get<2>(*it)+(scaling_step_delta/2), angle);
        if( corr2 > corr) {
          best_angle = angle;
          best_scale = std::get<2>(*it)+(scaling_step_delta/2);
          corr = corr2;
        }

        corr2 = Image::ColorImage::bc_invariant_correlation( main_image, templates[std::get<1>(*it)],
                                       i, std::get<0>(*it), std::get<2>(*it)-(scaling_step_delta/2), angle);
        if( corr2 > corr) {
          best_angle = angle;
          best_scale = std::get<2>(*it)-(scaling_step_delta/2);
          corr = corr2;
        }

        corr2 = Image::ColorImage::bc_invariant_correlation( main_image, templates[std::get<1>(*it)],
                                       i, std::get<0>(*it), std::get<2>(*it)-(scaling_step_delta/2), angle-(rotation_step_delta/2));
        if( corr2 > corr) {
          best_scale = std::get<2>(*it)-(scaling_step_delta/2);
          best_angle = angle-(rotation_step_delta/2);
          corr = corr2;
        }

        corr2 = Image::ColorImage::bc_invariant_correlation( main_image, templates[std::get<1>(*it)],
                                       i, std::get<0>(*it), std::get<2>(*it)-(scaling_step_delta/2), angle+(rotation_step_delta/2));
        if( corr2 > corr) {
          best_scale = std::get<2>(*it)-(scaling_step_delta/2);
          best_angle = angle+(rotation_step_delta/2);
          corr = corr2;
        }

        corr2 = Image::ColorImage::bc_invariant_correlation( main_image, templates[std::get<1>(*it)],
                                       i, std::get<0>(*it), std::get<2>(*it)+(scaling_step_delta/2), angle-(rotation_step_delta/2));
        if( corr2 > corr) {
          best_scale = std::get<2>(*it)+(scaling_step_delta/2);
          best_angle = angle-(rotation_step_delta/2);
          corr = corr2;
        }

        corr2 = Image::ColorImage::bc_invariant_correlation( main_image, templates[std::get<1>(*it)],
                                       i, std::get<0>(*it), std::get<2>(*it)+(scaling_step_delta/2), angle+(rotation_step_delta/2));
        if( corr2 > corr) {
          best_scale = std::get<2>(*it)+(scaling_step_delta/2);
          best_angle = angle+(rotation_step_delta/2);
          corr = corr2;
        }

        if( corr < th3)
          continue;

#if SHOW_FILTERS == 1
#pragma omp critical (third)
{
        third_grade_pixels.push_back( std::make_tuple( i, std::get<0>(*it)));
}
#endif

        int _id = templates[std::get<1>(*it)].get_id();
        unsigned int height = static_cast<unsigned int>(round(templates[std::get<1>(*it)].get_height() * best_scale));
        unsigned int width = static_cast<unsigned int>(round(templates[std::get<1>(*it)].get_width() * best_scale));
        pDsCell dscell = new DisjointSet::DsCell< std::tuple<int,int,int,fp,uint,uint,float> >(
                                  std::make_tuple(_id, static_cast<int>(std::get<0>(*it)), static_cast<int>(i), corr, height, width, best_angle));
#pragma omp critical (res)
{
        results.push_back( dscell);
}

      }

    } // i

    free(main_ras_l); free(main_ras_a); free(main_ras_b);
    free(buff_l);
    free(buff_l_S);
    free(buff_l_S2);
    free(aux); free(aux2);
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

  Image::ColorImage mask_image3( main_image);
  for( std::vector< std::tuple< unsigned int, unsigned int> >::iterator it = third_grade_pixels.begin();
       it != third_grade_pixels.end(); ++it) {
    mask_image3.L(std::get<0>(*it), std::get<1>(*it)) = LMAGENTA;
    mask_image3.A(std::get<0>(*it), std::get<1>(*it)) = AMAGENTA;
    mask_image3.B(std::get<0>(*it), std::get<1>(*it)) = BMAGENTA;
  }

  Image::ColorImage::write_image_to_bitmap( mask_image3, "m3_.bmp");
#endif

  /* SLINK clustering */
  unsigned int res_n = results.size();
  fp min_dist = pow( min_radius, 2);
  for(i=0; i<res_n; i++) {
    pDsCell ci = results[i];
    for(j=i+1; j<res_n; j++) {
      pDsCell cj = results[j];

      fp dist = pow( std::get<1>(ci->data) - std::get<1>(cj->data), 2) +
                pow( std::get<2>(ci->data) - std::get<2>(cj->data), 2);
      if( dist <= min_dist)
        DisjointSet::ds_union( *ci, *cj);
    }
  }

  /* move best correlation pixel to cluster parent */
  std::vector< pDsCell > best_results;
  for( std::vector< pDsCell >::iterator it = results.begin(); it != results.end(); ++it) {
    pDsCell root = &DisjointSet::ds_find(*(*it));
    if( (*it) == root ) {
      best_results.push_back(*it);
    }
    else if( std::get<3>((*it)->data) > std::get<3>(root->data) ) {
      std::swap( (*it)->data , root->data);
    }
  }

#if FRAME_TARGET==1
  Image::ColorImage target_main_image( main_image);
#endif

  /* print cluster parents only */
  for( std::vector< pDsCell >::iterator it = best_results.begin(); it != best_results.end(); ++it) {
    float angle = std::get<6>((*it)->data);
    int h = std::get<4>((*it)->data) / 2;
    int w = std::get<5>((*it)->data) / 2;
    int y = std::get<2>((*it)->data) - static_cast<int>(round(w*sin( Utils::D2R * angle) + h*cos( Utils::D2R * angle)));
    int x = std::get<1>((*it)->data) - static_cast<int>(round(w*cos( Utils::D2R * angle) - h*sin( Utils::D2R * angle)));
    std::cout << std::get<0>((*it)->data) << '\t' << x << '\t' << y << std::endl;
#if FRAME_TARGET==1
    Image::frame_target( y, x, h*2, w*2, angle, target_main_image);
#endif
  }

#if FRAME_TARGET==1
  Image::ColorImage::write_image_to_bitmap( target_main_image, "target.bmp");
#endif

  for( std::vector< pDsCell >::iterator it = results.begin(); it != results.end(); ++it)
    delete (*it);

  return 0;
}
