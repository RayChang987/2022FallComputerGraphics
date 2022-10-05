///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include<unordered_map>
#include <time.h>
using namespace std;

// constants
const int           RED = 0;                // red channel
const int           GREEN = 1;                // green channel
const int           BLUE = 2;                // blue channel
const unsigned char BACKGROUND[3] = { 0, 0, 0 };      // background color

// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1; i <= s; i++)
        res = (n - i + 1) * res / i;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
    data = new unsigned char[width * height * 4];
    ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char* d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
        data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image)
{
    width = image.width;
    height = image.height;
    data = NULL;
    if (image.data != NULL) {
        data = new unsigned char[width * height * 4];
        memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char* rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (!data)
        return NULL;

    // Divide out the alpha
    for (i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        int out_offset = i * width * 3;

        for (j = 0; j < width; j++)
        {
            RGBA_To_RGB(data + (in_offset + j * 4), rgb + (out_offset + j * 3));
        }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char* filename)
{
    TargaImage* out_image = Reverse_Rows();

    if (!out_image)
        return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
        cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
        return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char* filename)
{
    unsigned char* temp_data;
    TargaImage* temp_image;
    TargaImage* result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
        width = height = 0;
        return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////


bool TargaImage::To_Grayscale()  // passed
{
    data_to_imgRGBA();

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            float r = img_RGBA[0][i][j], g = img_RGBA[1][i][j], b = img_RGBA[2][i][j];
            unsigned char gray_scale_value = 0.299 * r + 0.587 * g + 0.114 * b;
            img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = gray_scale_value;
        }
    }
    imgRGBA_to_data();
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform() //tested
{
    data_to_imgRGBA();
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            img_RGBA[0][i][j] = ((img_RGBA[0][i][j] >> 5) << 5);
            img_RGBA[1][i][j] = ((img_RGBA[1][i][j] >> 5) << 5);
            img_RGBA[2][i][j] = ((img_RGBA[2][i][j] >> 6) << 6);
        }
    }
    imgRGBA_to_data();
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
int L2(int r1, int r2, int b1, int b2, int g1, int g2) {
    return (r2 - r1) * (r2 - r1) + (b2 - b1) * (b2 - b1) + (g2 - g1) * (g2 - g1);
}
bool TargaImage::Quant_Populosity() //tested
{
    data_to_imgRGBA();
    map<int, int> feq;
    using pixel = vector<int>;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int encoded_data = ((img_RGBA[0][i][j] >> 3) << 10) + ((img_RGBA[1][i][j] >> 3) << 5) + ((img_RGBA[2][i][j] >> 3) << 0);
            feq[encoded_data]++;
        }
    }
    vector<pair<int, int>> colors(feq.begin(), feq.end());
    vector<pixel> decoded_colors;
    sort(colors.begin(), colors.end(), [](const pair<int, int>& a, const  pair<int, int>& b) {return a.second > b.second;});
    int colorN = 256;
    for (int i = 0; i < colorN; ++i) {
        int r = (colors[i].first >> 10);
        colors[i].first -= (r << 10);
        int g = (colors[i].first >> 5);
        colors[i].first -= (g << 5);
        int b = colors[i].first;
        decoded_colors.push_back({ (r << 3), (g << 3), (b << 3) });
    }
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int mn = 1e9, mnI = -1;
            for (int k = 0; k < colorN; ++k) {
                int dis = L2(img_RGBA[0][i][j], decoded_colors[k][0], img_RGBA[1][i][j], decoded_colors[k][1], img_RGBA[2][i][j], decoded_colors[k][2]);
                if (dis < mn) {
                    mnI = k;
                    mn = dis;
                }
            }
            img_RGBA[0][i][j] = decoded_colors[mnI][0];
            img_RGBA[1][i][j] = decoded_colors[mnI][1];
            img_RGBA[2][i][j] = decoded_colors[mnI][2];
        }
    }
    //for (auto& v : colors) {
    //    cout << (int)v[0] << " " << (int)v[1] << " " << (int)v[2] << endl;
    //}
    imgRGBA_to_data();
    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()//tested
{
    data_to_imgRGBA();
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            float r = img_RGBA[0][i][j], g = img_RGBA[1][i][j], b = img_RGBA[2][i][j];
            float gray_scale_value = floor(0.299 * r + 0.587 * g + 0.114 * b) / 256;
            if ((gray_scale_value) >= 0.5) {
                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = 255;
            }
            else {
                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = 0;
            }

        }
    }
    imgRGBA_to_data();
    return true;
}// Dither_Threshold

bool TargaImage::Compare(TargaImage* pImage) {
    bool re = true;
    if (width != pImage->width) {
        return false;
    }
    if (height != pImage->height) {
        return false;
    }
    for (int i = 0; i < width * height * 4; i++) {
        //if (i % 3 == 0) { continue; }
        if (data[i] != pImage->data[i]) {
            //return false;
            re = false;
            data[i - (i & 3)] = 255;
            data[i - (i & 3) + 1] = 0;
            data[i - (i & 3) + 2] = 0;
        }
    }
    return re;
}
///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()//tested
{
    srand(time(NULL));

    data_to_imgRGBA();
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            float r = img_RGBA[0][i][j], g = img_RGBA[1][i][j], b = img_RGBA[2][i][j];
            float gray_scale_value = floor(0.299 * r + 0.587 * g + 0.114 * b) / 256;
            int rand_val = rand() % 1000;
            float x = (float)((0.4) * rand_val) / 1000 - 0.2;
            gray_scale_value += x;
            if (gray_scale_value >= 0.5) {
                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = 255;
            }
            else {
                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = 0;
            }

        }
    }
    imgRGBA_to_data();
    return true;

}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS() //tested
{
    data_to_imgRGBA();
    vector<vector<float>> d(height, vector<float>(width));
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            float r = img_RGBA[0][i][j], g = img_RGBA[1][i][j], b = img_RGBA[2][i][j];
            unsigned char gray_scale_value = 0.299 * r + 0.587 * g + 0.114 * b;
            img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = gray_scale_value;
            d[i][j] = gray_scale_value;
        }
    }
    
    for (int i = 0; i < height; i++){
        if (i % 2 == 0){
            for (int j = 0; j < width; j++)
            {
                float result = (d[i][j] >= 127.5) ? 255 : 0;
                float error = d[i][j] - result;

                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = result;
                img_RGBA[3][i][j] = 255;
                if (i + 1 < height && j > 0) {
                    d[i + 1][j - 1] += error * ((float)3 / 16);
                }
                if (i + 1 < height) {
                    d[i + 1][j] += error * ((float)5 / 16);
                }
                if (i + 1 < height && j + 1 < width) {
                    d[i + 1][j + 1] += error * ((float)1 / 16);
                }
                if (j + 1 < width) {
                    d[i][j + 1] += error * ((float)7 / 16);
                }
                
            }
        }
        else{
            for (int j = width - 1; j >= 0; j--){
                float result = (d[i][j] >= 127.5) ? 255 : 0;
                float error = d[i][j] - result;

                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = result;
                img_RGBA[3][i][j] = 255;
                if (i+1 < height && j+1 < width) {
                    d[i + 1] [j + 1] += error * ((float)3 / 16);
                }
                if (i+1 < height) {
                    d[i + 1] [j] += error * ((float)5 / 16);
                }
                if (i+1< height && j>0) {
                    d[i + 1] [j - 1] += error * ((float)1 / 16);
                }
                if (j > 0) {
                    d[i] [j - 1] += error * ((float)7 / 16);
                }
            }
        }
    }
    imgRGBA_to_data();
    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()//tested
{
    data_to_imgRGBA();
    float sum = 0;
    vector<float> gray_scale_values = { 255 };
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            float r = img_RGBA[0][i][j], g = img_RGBA[1][i][j], b = img_RGBA[2][i][j];
            int gray_scale_value = (0.299 * r + 0.587 * g + 0.114 * b);
            sum += gray_scale_value;
            gray_scale_values.push_back(gray_scale_value);
        }
    }
    //cout << "The average Brightness is: " << (sum / (height * width)) << endl;
    sort(gray_scale_values.begin(), gray_scale_values.end());
    reverse(gray_scale_values.begin(), gray_scale_values.end());
    float average = sum / (float)(height * width);

    int number_of_1 = (int)(sum / 256);
    float threshold = gray_scale_values[number_of_1] / 256;
    cout << "Using Threshold : " << threshold << endl;
    float cnt = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            float r = img_RGBA[0][i][j], g = img_RGBA[1][i][j], b = img_RGBA[2][i][j];
            float gray_scale_value = floor(0.299 * r + 0.587 * g + 0.114 * b) / 256;
            if (gray_scale_value >= threshold) {
                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = 255;
                cnt += 255;
            }
            else {
                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = 0;
            }
        }
    }
    //0.176854
    imgRGBA_to_data();
    cout << "The average Brightness before is: " << (sum / (height * width)) << "-->" << (sum / (height * width) / 256) << endl;
    cout << "The average Brightness after is: " << (cnt / (height * width)) << "-->" << (cnt / (height * width) / 256) << endl;
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()//tested
{
    data_to_imgRGBA();
    float mask[4][4] = {
                    {0.7059,0.3529,0.5882,0.2353},
                    {0.0588,0.9412,0.8235,0.4118},
                    {0.4706,0.7647,0.8824,0.1176},
                    {0.1765,0.5294,0.2941,0.6471}
    };
   
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            float r = img_RGBA[0][i][j], g = img_RGBA[1][i][j], b = img_RGBA[2][i][j];
            float gray_scale_value = (0.299 * r + 0.587 * g + 0.114 * b) / 256.0;
            float msk = mask[(i) % 4][(j) % 4];
            if (gray_scale_value >= msk) {
                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = 255;
            }
            else {
                img_RGBA[0][i][j] = img_RGBA[1][i][j] = img_RGBA[2][i][j] = 0;
            }
        }
    }
    imgRGBA_to_data();
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
vector<int> _3_bit_color = { 0, 36, 73, 109, 146, 182, 219, 235 };
vector<int> _2_bit_color = { 0, 85, 170, 255 };
int _3_bit_color_closest(int x) {
    int mn = INT_MAX, mnI = -1;
    for (int i = 0; i < 8; ++i) {
        int d = abs(_3_bit_color[i] - x);
        if (d < mn) {
            mn = d;
            mnI = i;
        }
    }
    if (mnI == -1) { cerr << x << endl; system("PAUSE"); }
    return mnI;
}
int _2_bit_color_closest(int x) {
    int mn = INT_MAX, mnI = -1;
    for (int i = 0; i < 4; ++i) {
        int d = abs(_2_bit_color[i] - x);
        if (d < mn) {
            mn = d;
            mnI = i;
        }
    }
    if (mnI == -1) { cerr << x << endl; system("PAUSE"); }
    return mnI;
}
bool TargaImage::Dither_Color()
{
    data_to_imgRGBA();
    vector<vector<vector<float>>> d(4, vector<vector<float>>(height, vector<float>(width)));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < height; ++j) {
            for (int k = 0; k < width; ++k) {
                d[i][j][k] = img_RGBA[i][j][k];
            }
        }
    }

    for (int k = 0; k < 3; ++k) {
        for (int i = 0; i < height; i++) {
            if (i % 2 == 0) {
                for (int j = 0; j < width; j++)
                {
                    int idx = (k == 2) ? _2_bit_color_closest(d[k][i][j]): _3_bit_color_closest(d[k][i][j]);
                    float cmp = (k == 2) ? _2_bit_color[idx] : _3_bit_color[idx];
                    float result;
                    if (k == 2) {
                        if (d[k][i][j]>=cmp) {result = cmp;}
                        else {result = _2_bit_color[idx - 1];}
                    }
                    else {
                        if (d[k][i][j] >= cmp) { result = cmp; }
                        else { result = _3_bit_color[idx - 1]; }
                    }
                    float error = d[k][i][j] - result;

                    img_RGBA[k][i][j] = result;
                    if (i + 1 < height && j > 0) {
                        d[k][i + 1][j - 1] += error * ((float)3 / 16);
                    }
                    if (i + 1 < height) {
                        d[k][i + 1][j] += error * ((float)5 / 16);
                    }
                    if (i + 1 < height && j + 1 < width) {
                        d[k][i + 1][j + 1] += error * ((float)1 / 16);
                    }
                    if (j + 1 < width) {
                        d[k][i][j + 1] += error * ((float)7 / 16);
                    }

                }
            }
            else {
                for (int j = width - 1; j >= 0; j--) {
                    int idx = (k == 2) ? _2_bit_color_closest(d[k][i][j]) : _3_bit_color_closest(d[k][i][j]);
                    float cmp = (k == 2) ? _2_bit_color[idx] : _3_bit_color[idx];
                    float result;
                    if (k == 2) {
                        if (d[k][i][j] >= cmp) { result = cmp; }
                        else { result = _2_bit_color[idx - 1]; }
                    }
                    else {
                        if (d[k][i][j] >= cmp) { result = cmp; }
                        else { result = _3_bit_color[idx - 1]; }
                    }
                    float error = d[k][i][j] - result;
                    if (i + 1 < height && j + 1 < width) {
                        d[k][i + 1][j + 1] += error * ((float)3 / 16);
                    }
                    if (i + 1 < height) {
                        d[k][i + 1][j] += error * ((float)5 / 16);
                    }
                    if (i + 1 < height && j > 0) {
                        d[k][i + 1][j - 1] += error * ((float)1 / 16);
                    }
                    if (j > 0) {
                        d[k][i][j - 1] += error * ((float)7 / 16);
                    }
                }
            }
        }
    }
    
    imgRGBA_to_data();
    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Over: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        cout << width << " " << pImage->width << " " << endl;
        return false;
    }// if

    for (int i = 0; i < width * height * 4; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i + 1] = abs(rgb1[1] - rgb2[1]);
        data[i + 2] = abs(rgb1[2] - rgb2[2]);
        data[i + 3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    data_to_imgRGBA();
    vector<vector<vector<int>>> t(4, vector<vector<int>>(height, vector<int>(width)));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < 4; ++k) {
                t[k][i][j] = img_RGBA[k][i][j];
            }
        }
    }
    float kernel[5][5] = {
        {1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1}
    };
    float d = 0;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            d += kernel[i][j];
        }
    }
   
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < 4; k++) {
                float sum = 0;
                for(int p = -2; p<=2; ++p){
                    for(int q = -2; q<=2; ++q){
                        int ni = i + p, nj = j + q;
                        if (ni < 0 || nj < 0 || ni >= height || nj >= width) {continue;}
                        sum += t[k][ni][nj] * kernel[p + 2][q + 2];
                    }
                }
                img_RGBA[k][i][j] = sum / d;
            }

        }
    }
    imgRGBA_to_data();
    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    data_to_imgRGBA();
    vector<vector<vector<int>>> t(4, vector<vector<int>>(height, vector<int>(width)));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < 4; ++k) {
                t[k][i][j] = img_RGBA[k][i][j];
            }
        }
    }
    float kernel[5][5] = {
        {1, 2, 3, 2, 1},
        {2, 4, 6, 4, 2},
        {3, 6, 9, 6, 3},
        {2, 4, 6, 4, 2},
        {1, 2, 3, 2, 1}
    };
    float d = 0;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            d += kernel[i][j];
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < 4; k++) {
                float sum = 0;
                for (int p = -2; p <= 2; ++p) {
                    for (int q = -2; q <= 2; ++q) {
                        int ni = i + p, nj = j + q;
                        if (ni < 0 || nj < 0 || ni >= height || nj >= width) { continue; }
                        sum += t[k][ni][nj] * kernel[p + 2][q + 2];
                    }
                }
                img_RGBA[k][i][j] = sum / d;
            }

        }
    }
    imgRGBA_to_data();
    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    data_to_imgRGBA();
    vector<vector<vector<int>>> t(4, vector<vector<int>>(height, vector<int>(width)));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < 4; ++k) {
                t[k][i][j] = img_RGBA[k][i][j];
            }
        }
    }
    float kernel[5][5] = {
        {1, 4, 6, 4, 1},
        {4, 16, 24, 16, 4},
        {6, 24, 36, 24, 6},
        {4, 16, 24, 16, 4},
        {1, 4, 6, 4, 1}
    };
    float d = 0;
   
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            d += kernel[i][j];
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < 4; k++) {
                float sum = 0;
                for (int p = -2; p <= 2; ++p) {
                    for (int q = -2; q <= 2; ++q) {
                        int ni = i + p, nj = j + q;
                        if (ni < 0 || nj < 0 || ni >= height || nj >= width) { continue; }
                        sum += t[k][ni][nj] * kernel[p + 2][q + 2];
                    }
                }
                img_RGBA[k][i][j] = sum / d;
            }

        }
    }
    imgRGBA_to_data();
    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N(unsigned int N)
{
    data_to_imgRGBA();
    vector<vector<vector<int>>> t(4, vector<vector<int>>(height, vector<int>(width)));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < 4; ++k) {
                t[k][i][j] = img_RGBA[k][i][j];
            }
        }
    }
    const int kernelN = N;
    vector<vector<float>> kernel(kernelN, vector<float>(kernelN));
    vector<int> v(kernelN);
    
    for (int i = 0; i < kernelN; ++i) {v[i] = Binomial(N-1, i);}
    for (int i = 0; i < kernelN; ++i) { cout << v[i] << " "; }
    for (int i = 0; i < kernelN; ++i) {
        for (int j = 0; j < kernelN; ++j) {
            kernel[i][j] = v[i] * v[j];
        }
    }
    cout << "Using Kernel : " << endl;
    for (int i = 0; i < kernelN; ++i) {
        for (int j = 0; j < kernelN; ++j) {
            cout << kernel[i][j] << " ";
        }
        cout << endl;
    }
    float d = 0;
    for (int i = 0; i < kernelN; ++i) {
        for (int j = 0; j < kernelN; ++j) {
            d += kernel[i][j];
        }
    }
    int r = kernelN / 2;
    for (int i = r; i < height-r; i++) {
        for (int j = r; j < width-r; j++) {
            for (int k = 0; k < 4; k++) {
                float sum = 0;
                for (int p = -r; p <= r; ++p) {
                    for (int q = -r; q <= r; ++q) {
                        int ni = i + p, nj = j + q;
                        if (ni < 0 || nj < 0 || ni >= height || nj >= width) { continue; } //not neccessary
                        sum += t[k][ni][nj] * kernel[p + r][q + r];
                    }
                }
                img_RGBA[k][i][j] = sum / d;
            }

        }
    }
    imgRGBA_to_data();
    return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    ClearToBlack();
    return false;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    ClearToBlack();
    return false;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    data_to_imgRGBA();
    vector<vector<vector<int>>> t(4, vector<vector<int>>(height, vector<int>(width)));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < 4; ++k) {
                t[k][i][j] = img_RGBA[k][i][j];
            }
        }
    }
    float kernel1[3][3] = {
        {1.0/16.0, 1.0/8.0, 1.0/16.0},
        {1.0/8.0, 1.0/4.0, 1.0/8.0}, 
        {1.0/16.0, 1.0/8.0, 1.0/16.0}
    };
    float d = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            d += kernel1[i][j];
        }
    }
    int new_height = height / 2, new_width = width / 2;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (i & 1 || j & 1) { continue; }
            for (int k = 0; k < 4; k++) {
                float sum = 0;
                for (int p = -1; p <= 1; ++p) {
                    for (int q = -1; q <= 1; ++q) {
                        int ni = i + p, nj = j + q;
                        if (ni < 0 || nj < 0 || ni >= height || nj >= width) { continue; }
                        sum += t[k][ni][nj] * kernel1[p + 1][q + 1];
                    }
                }
                img_RGBA[k][i/2][j/2] = sum / d;
            }

        }
    }
    height = new_height;
    width = new_width;
    imgRGBA_to_data();
    return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    data_to_imgRGBA();
    vector<vector<vector<int>>> t(4, vector<vector<int>>(height, vector<int>(width)));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < 4; ++k) {
                t[k][i][j] = img_RGBA[k][i][j];
            }
        }
    }
    float kernel1[3][3] = {
        {1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0},
        {1.0 / 8.0, 1.0 / 4.0, 1.0 / 8.0},
        {1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0}
    };
    float kernel2[4][4] = {
        {1.0 / 64.0, 3.0 / 64.0, 3.0 / 64.0, 1.0 / 64.0},
        {3.0 / 64.0, 9.0 / 64.0, 9.0 / 64.0, 3.0 / 64.0},
        {3.0 / 64.0, 9.0 / 64.0, 9.0 / 64.0, 3.0 / 64.0},
        {1.0 / 64.0, 3.0 / 64.0, 3.0 / 64.0, 1.0 / 64.0}
    };
    float kernel3[4][3] = {
        { 1.0 / 32.0, 2.0 / 32.0, 1.0/ 32.0},
        { 3.0 / 32.0, 6.0/ 32.0, 3.0 / 32.0},
        { 3.0 / 32.0, 6.0 / 32.0, 3.0 / 32.0},
        { 1.0 / 32.0, 2.0 / 32.0, 1.0 / 32.0}
    };
    float d1 = 0, d2 = 0, d3 = 0;
    int kN1 = 3, kN2 = 4, r1 = kN1/2, r2 = kN2/2;
    for (int i = 0; i < kN1; ++i) {
        for (int j = 0; j < kN1; ++j) {
            d1 += kernel1[i][j];
        }
    }
    for (int i = 0; i < kN2; ++i) {
        for(int j = 0; j<kN2; ++j){
            d2+=kernel2[i][j];
        }
    }
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            d3 += kernel3[i][j];
        }
    }
    int new_height = height*2, new_width = width*2;
    //cout << r1 << " " << r2 << endl;
    //cout << d1 << " " << d2 << " " << d3 << endl;
    for (auto& mat : img_RGBA) {
        mat.resize(new_height);
        for (auto& row : mat) {
            row.resize(new_width);
        }
    }
    
    for (int i = 0; i < new_height; i++) {
        for (int j = 0; j < new_width; j++) {
            for (int k = 0; k < 4; k++) {
                float sum = 0;
                if (i % 2 == 0 && j % 2 == 0) {
                    for (int p = -1; p <= 1; ++p) {
                        for (int q = -1; q <= 1; ++q) {
                            int ni = i + p, nj = j + q;
                            if (ni < 0 || nj < 0 || ni >= new_height || nj >= new_width) { continue; }
                            sum += t[k][ni / 2][nj / 2] * kernel1[p + r1][q + r1];
                        }
                    }
                    img_RGBA[k][i][j] = sum / d1;
                }
                else  if (i%2==1&&j%2==1) {
                    for (int p = -1; p <= 2; ++p) {
                        for (int q = -1; q <= 2; ++q) {
                            int ni = i + p, nj = j + q;
                            if (ni < 0 || nj < 0 || ni >= new_height || nj >= new_width) { continue; }
                            sum += t[k][ni / 2][nj / 2] * kernel2[p + r1][q + r1];
                        }
                    }
                    img_RGBA[k][i][j] = sum / d2;
                }
                else if(i%2==1&&j%2==0){
                    for (int p = -2; p <= 1; ++p) {
                        for (int q = -1; q <= 1; ++q) {
                            int ni = i + p, nj = j + q;
                            if (ni < 0 || nj < 0 || ni >= new_height || nj >= new_width) { continue; }
                            sum += t[k][ni / 2][nj / 2] * kernel3[p +2][q + 1];
                        }
                    }
                    img_RGBA[k][i][j] = sum / d3;
                }
                else if (i % 2 == 0 && j % 2 == 1) {
                    for (int p = -2; p <= 1; ++p) {
                        for (int q = -1; q <= 1; ++q) {
                            int ni = i + q, nj = j + p;
                            if (ni < 0 || nj < 0 || ni >= new_height || nj >= new_width) { continue; }
                            sum += t[k][ni / 2][nj / 2] * kernel3[p + 2][q + 1];
                        }
                    }
                    img_RGBA[k][i][j] = sum / d3;
                }
                //for (int p = -1; p <= 1; ++p) {
                //    for (int q = -1; q <= 1; ++q) {
                //        int ni = i + p, nj = j + q;
                //        if (ni < 0 || nj < 0 || ni >= new_height || nj >= new_width) { continue; }
                //        sum += t[k][ni / 2][nj / 2] * kernel1[p + r1][q + r1];
                //    }
                //}
                //img_RGBA[k][i][j] = sum / d1;

            }

        }
    }
    
    
    unsigned char* new_data = new unsigned char[width*2 * height*2 * 4];
    for (int i = 0; i < height * width * 4; ++i) {
        new_data[i] = data[i];
    }
    height = new_height; width = new_width;
    swap(new_data, data);
    delete[]new_data;
    //TODO::delete new_data
    new_data = NULL;
    
    imgRGBA_to_data();
    return true;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    data_to_imgRGBA();
    vector<vector<vector<int>>> t(4, vector<vector<int>>(height, vector<int>(width)));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < 4; ++k) {
                t[k][i][j] = img_RGBA[k][i][j];
            }
        }
    }
    float kernel1[3][3] = {
        {1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0},
        {1.0 / 8.0, 1.0 / 4.0, 1.0 / 8.0},
        {1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0}
    };

    imgRGBA_to_data();
    return true;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Rotate(float angleDegrees)
{
    const float PI = 3.1415926;
    angleDegrees = angleDegrees * PI / 180.0;
    data_to_imgRGBA();
    vector<vector<vector<int>>> t(4, vector<vector<int>>(height, vector<int>(width)));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < 4; ++k) {
                t[k][i][j] = img_RGBA[k][i][j];
            }
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int ii = int(cos(angleDegrees) * (i- (height/2)) - sin(angleDegrees) * (j-(width/2)));
            int jj = int(sin(angleDegrees) * (i - (height/2)) + cos(angleDegrees) * (j-(width/2)));
            ii += height / 2, jj += width / 2;
            //cout << rotated_i << " " << rotated_j << endl;
            if (ii < 0 || ii >= height || jj < 0 || jj >= width) {
                for (int k = 0; k < 4; k++) { img_RGBA[k][i][j] = 0; }
            }
            else{
                for (int k = 0; k < 4; k++) { img_RGBA[k][i][j] = t[k][ii][jj]; }
            }
            

        }
    }
    vector<vector<vector<int>>> t2(4, vector<vector<int>>(height, vector<int>(width)));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < 4; ++k) {
                t2[k][i][j] = img_RGBA[k][i][j];
            }
        }
    }
	float kernel[3][3] = {
		{1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0},
		{1.0 / 8.0, 1.0 / 4.0, 1.0 / 8.0},
		{1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0}
	};
    float d = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            d += kernel[i][j];
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < 4; k++) {
                float sum = 0;
                for (int p = -1; p < 1; ++p) {
                    for (int q = -1; q < 1; ++q) {
                        int ni = i + p, nj = j + q;
                        if (ni < 0 || nj < 0 || ni >= height || nj >= width) { continue; }
                        sum += t2[k][ni][nj] * kernel[p + 1][q + 1];
                    }
                }
                img_RGBA[k][i][j] = sum / d;
            }

        }
    }
    imgRGBA_to_data();
    return true;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//`
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char* rgba, unsigned char* rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
        float	alpha_scale = (float)255 / (float)alpha;
        int	val;
        int	i;

        for (i = 0; i < 3; i++)
        {
            val = (int)floor(rgba[i] * alpha_scale);
            if (val < 0)
                rgb[i] = 0;
            else if (val > 255)
                rgb[i] = 255;
            else
                rgb[i] = val;
        }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char* dest = new unsigned char[width * height * 4];
    TargaImage* result;
    int 	        i, j;

    if (!data)
        return NULL;

    for (i = 0; i < height; i++)
    {
        int in_offset = (height - i - 1) * width * 4;
        int out_offset = i * width * 4;

        for (j = 0; j < width; j++)
        {
            dest[out_offset + j * 4] = data[in_offset + j * 4];
            dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
            dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
            dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
    int radius_squared = (int)s.radius * (int)s.radius;
    for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
        for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
            int x_loc = (int)s.x + x_off;
            int y_loc = (int)s.y + y_off;
            // are we inside the circle, and inside the image?
            if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
                int dist_squared = x_off * x_off + y_off * y_off;
                if (dist_squared <= radius_squared) {
                    data[(y_loc * width + x_loc) * 4 + 0] = s.r;
                    data[(y_loc * width + x_loc) * 4 + 1] = s.g;
                    data[(y_loc * width + x_loc) * 4 + 2] = s.b;
                    data[(y_loc * width + x_loc) * 4 + 3] = s.a;
                }
                else if (dist_squared == radius_squared + 1) {
                    data[(y_loc * width + x_loc) * 4 + 0] =
                        (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
                    data[(y_loc * width + x_loc) * 4 + 1] =
                        (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
                    data[(y_loc * width + x_loc) * 4 + 2] =
                        (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
                    data[(y_loc * width + x_loc) * 4 + 3] =
                        (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
    unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
    radius(iradius), x(ix), y(iy), r(ir), g(ig), b(ib), a(ia)
{
}

void TargaImage::data_to_imgRGBA() {
    img_RGBA.resize(4); //R, G, B, A
    for (auto& m : img_RGBA) {
        m.resize(height);
        for (auto& row : m) {
            row.resize(width);
        }
    }

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int offset = i * width * 4 + j * 4;
            for (int k = 0; k < 4; ++k) {
                img_RGBA[k][i][j] = data[offset + k];
            }
        }
    }
}
void TargaImage::imgRGBA_to_data() {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int  offset = i * width * 4 + j * 4;
            for (int k = 0; k < 4; ++k) {
                data[offset + k] = img_RGBA[k][i][j];
            }
        }
    }
}


mat TargaImage::mul(mat& a, mat& b) {
    if (a[0].size() != b.size()) { throw("Shape Error!"); }
    int m = a.size(), n = a[0].size(), p = b[0].size();
    mat ret(m, vector<unsigned char>(p));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            ret[i][j] = 0;
            for (int k = 0; k < n; k++) {
                ret[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return ret;
}
