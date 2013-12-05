#include <iostream>
#include "colorimage.h"

using namespace std;

int main()
{
    Image::ColorImage img("test_case_2/4sequoia_test.bmp");
    Image::ColorImage::write_image_to_bitmap(img.scale_image(1.2f), "scaled12.bmp");
    Image::ColorImage::write_image_to_bitmap(img.scale_image(0.5f), "scaled05.bmp");
    Image::ColorImage temp("test_case_2/001template.bmp");
    Image::ColorImage::tag( img, temp, 500, 500, 2.f, 0.f);
    fp corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 2.f, 0.f);
    cout << corr << endl;

    cout<<endl;

    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.9f, 0.f);
    cout << corr << endl;
    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.8f, 0.f);
    cout << corr << endl;
    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.7f, 0.f);
    cout << corr << endl;

    cout<<endl;

    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 498, 499, 2.f, 0.f);
    cout << corr << endl;

    cout<<endl;

    corr = Image::ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.7f, 15.f);
    cout << corr << endl;

    Image::ColorImage::write_image_to_bitmap(img, "tagged.bmp");
    return 0;
}

//    Image::ColorImage::write_image_to_bitmap(img, "/home/albert/workspace/ayc/original.bmp");
//    Image::ColorImage blur = Image::ColorImage::gaussian_smoother( img);
//    Image::ColorImage::write_image_to_bitmap(blur, "/home/albert/workspace/ayc/blur.bmp");
//    Image::ColorImage::write_image_to_bitmap(img.scale_image(1.2f), "/home/albert/workspace/ayc/scaled12.bmp");
//    Image::ColorImage::write_image_to_bitmap(img.scale_image(0.5f), "/home/albert/workspace/ayc/scaled05.bmp");
//    Image::ColorImage::write_image_to_bitmap(img.rotate_image(10.f), "/home/albert/workspace/ayc/rotate10.bmp");
//    Image::ColorImage::write_image_to_bitmap(img.rotate_image(20.f), "/home/albert/workspace/ayc/rotate20.bmp");
