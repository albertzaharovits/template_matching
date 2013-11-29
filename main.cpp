#include <iostream>
#include "colorimage.h"

using namespace std;

int main()
{
    ColorImage img("test_case_2/4sequoia_test.bmp");
    ColorImage::write_image_to_bitmap(img.scale_image(1.2f), "scaled12.bmp");
    ColorImage::write_image_to_bitmap(img.scale_image(0.5f), "scaled05.bmp");
    ColorImage temp("test_case_2/001template.bmp");
    ColorImage::tag( img, temp, 500, 500, 2.f, 0.f);
    fp corr = ColorImage::bc_invariant_correlation( img, temp, 500, 500, 2.f, 0.f);
    cout << corr << endl;

    cout<<endl;

    corr = ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.9f, 0.f);
    cout << corr << endl;
    corr = ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.8f, 0.f);
    cout << corr << endl;
    corr = ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.7f, 0.f);
    cout << corr << endl;

    cout<<endl;

    corr = ColorImage::bc_invariant_correlation( img, temp, 498, 499, 2.f, 0.f);
    cout << corr << endl;

    cout<<endl;

    corr = ColorImage::bc_invariant_correlation( img, temp, 500, 500, 1.7f, 15.f);
    cout << corr << endl;

    ColorImage::write_image_to_bitmap(img, "tagged.bmp");
    return 0;
}

//    ColorImage::write_image_to_bitmap(img, "/home/albert/workspace/ayc/original.bmp");
//    ColorImage blur = ColorImage::gaussian_smoother( img);
//    ColorImage::write_image_to_bitmap(blur, "/home/albert/workspace/ayc/blur.bmp");
//    ColorImage::write_image_to_bitmap(img.scale_image(1.2f), "/home/albert/workspace/ayc/scaled12.bmp");
//    ColorImage::write_image_to_bitmap(img.scale_image(0.5f), "/home/albert/workspace/ayc/scaled05.bmp");
//    ColorImage::write_image_to_bitmap(img.rotate_image(10.f), "/home/albert/workspace/ayc/rotate10.bmp");
//    ColorImage::write_image_to_bitmap(img.rotate_image(20.f), "/home/albert/workspace/ayc/rotate20.bmp");
