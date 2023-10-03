///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
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
bool TargaImage::To_Grayscale()
{
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	for (int i = 0; i < (width * height * 4); i += 4) {

		unsigned char rgb[3];
		RGBA_To_RGB(data + i, rgb);
		unsigned char gray = (unsigned char)(0.299 * (float)rgb[0] + 0.587 * (float)rgb[1] + 0.114 * (float)rgb[2]);
		data[i + 0] = data[i + 1] = data[i + 2] = gray;
	}

	return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}

	for (int i = 0; i < (width * height * 4); i += 4) {
		unsigned char rgb[3];
		RGBA_To_RGB(data + i, rgb);
		//Reduce to 8 bits
		data[i + 0] = rgb[0] & (~((1 << (8 - 3)) - 1));//8 levels red
		data[i + 1] = rgb[1] & (~((1 << (8 - 3)) - 1));//8 levels green
		data[i + 2] = rgb[2] & (~((1 << (8 - 2)) - 1));//4 levels blue

	}
	//ClearToBlack();
	return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
	ClearToBlack();
	return false;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	To_Grayscale();
	for (int i = 0; i < width * height * 4; i += 4) {
		unsigned char rgb[3];
		RGBA_To_RGB(data + i, rgb);
		unsigned char result;
		if (rgb[0] > 128) {
			result = 255;
		}
		else {
			result = 0;
		}
		data[i + 0] = result;
		data[i + 1] = result;
		data[i + 2] = result;
	}
	return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}

	To_Grayscale();
	srand((unsigned int)time(NULL));

	float threashold = 0.5;
	for (int i = 0; i < width * height * 4; i += 4) {
		unsigned char rgb[3];
		RGBA_To_RGB(data + i, rgb);
		float random = ((float(rand()) / float(RAND_MAX)) * (float)(0.2 - (-0.2)) + (float)(-0.2));
		float result = (rgb[0] / (float)256) + random;
		if (result < threashold) {
			data[i + 0] = data[i + 1] = data[i + 2] = 0;
		}
		else {
			data[i + 0] = data[i + 1] = data[i + 2] = 255;
		}
	}
	return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
	ClearToBlack();
	return false;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	To_Grayscale();
	vector<unsigned char> image;
	float total = 0;
	for (int i = 0; i < width * height * 4; i += 4) {
		unsigned char rgb[3];
		RGBA_To_RGB(data + i, rgb);
		total += rgb[0];
		image.push_back(rgb[0]);
	}
	float average = (total / (width * height) / (float)255);
	float n = (1 - average) * (width * height);
	sort(image.begin(), image.end());
	unsigned char threshold = image[(unsigned int)n];
	for (int i = 0; i < width * height * 4; i += 4) {
		unsigned char rgb[3];
		RGBA_To_RGB(data + i, rgb);
		if (rgb[0] < threshold) {
			data[i + 0] = data[i + 1] = data[i + 2] = 0;
		}
		else {
			data[i + 0] = data[i + 1] = data[i + 2] = 255;
		}
	}
	return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	To_Grayscale();
	float matrix[4][4] = { {.75,    .375,   .625,   .25},{.0625,     1,   .875, .4375},{.5,    .8125,  .9375, .1250},{.1875, .5625,  .3125, .6875} };
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			unsigned char rgb[3];
			RGBA_To_RGB(data + (((x * width) + y) * 4), rgb);
			float result = (rgb[0] / (float)256);
			if (result >= matrix[x % 4][y % 4]) {
				data[(((x * width) + y) * 4) + 0] = data[(((x * width) + y) * 4) + 1] = data[(((x * width) + y) * 4) + 2] = 255;
			}
			else {
				data[(((x * width) + y) * 4) + 0] = data[(((x * width) + y) * 4) + 1] = data[(((x * width) + y) * 4) + 2] = 0;
			}
		}
	}
	return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
	ClearToBlack();
	return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
	if (!pImage)
		return false;
	if (width != pImage->width || height != pImage->height) {
		return false;
	}
	for (int i = 0; i < width * height * 4; i += 4) {
		double alpha = ((double)data[i + 3]) / 255.0;
		for (int j = 0; j < 4; j++)
			data[i + j] += (pImage->data[i + j] * (1.0 - alpha));
	}
	return true;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
	if (!pImage)
		return false;
	if (width != pImage->width || height != pImage->height) {
		return false;
	}
	for (int i = 0; i < width * height * 4; i += 4) {
		double alpha = ((double)pImage->data[i + 3] / 255.0);
		for (int j = 0; j < 4; j++) {
			data[i + j] = data[i + j] * alpha;
		}
	}
	return true;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
	if (!pImage)
		return false;
	if (width != pImage->width || height != pImage->height) {
		return false;
	}
	for (int i = 0; i < width * height * 4; i += 4) {
		double alpha = ((double)pImage->data[i + 3] / 255.0);
		for (int j = 0; j < 4; j++) {
			data[i + j] = data[i + j] * (1.0 - alpha);
		}
	}
	return true;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
	if (!pImage)
		return false;
	if (width != pImage->width || height != pImage->height) {
		return false;
	}
	for (int i = 0; i < width * height * 4; i += 4) {
		double aOne = ((double)data[i + 3] / 255.0);
		double aTwo = ((double)pImage->data[i + 3] / 255.0);
		for (int j = 0; j < 4; j++) {
			data[i + j] = (data[i + j] * aTwo) + (pImage->data[i + j] * (1.0 - aOne));
		}
	}
	return true;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
	if (!pImage)
		return false;
	if (width != pImage->width || height != pImage->height) {
		return false;
	}
	for (int i = 0; i < width * height * 4; i += 4) {
		double aOne = ((double)data[i + 3] / 255.0);
		double aTwo = ((double)pImage->data[i + 3] / 255.0);
		for (int j = 0; j < 4; j++) {
			data[i + j] = (data[i + j] * (1.0 - aTwo)) + (pImage->data[i + j] * (1.0 - aOne));
		}
	}
	return true;
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

	if (width != pImage->width || height != pImage->height) {
		cout << "Difference: Images not the same size\n";
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
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	double box[5][5];
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			box[i][j] = 1.0 / 25.0;
		}
	}
	unsigned char* rgb = To_RGB();
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width; y++) {
			for (int color = 0; color < 3; color++) {
				double total = 0;
				for (int r = 0; r < 5; r++) {
					for (int c = 0; c < 5; c++) {
						int rPos = x - 2 + r;
						int cPos = y - 2 + c;
						if (rPos < 0) {
							rPos = -rPos;
						}
						if (rPos >= height) {
							rPos = ((height * 2) - 2) - rPos;
						}
						if (cPos < 0) {
							cPos = -cPos;
						}
						if (cPos >= width) {
							cPos = ((width * 2) - 2) - cPos;
						}
						total += rgb[((rPos * width + cPos) * 3) + color] * box[r][c];
					}
				}
				data[(x * width + y) * 4 + color] = (unsigned char)total;
				if (total > 255) {
					data[(x * width + y) * 4 + color] = 255;
				}
				else if (total < 0) {
					data[(x * width + y) * 4 + color] = 0;
				}
			}
		}
	}
	delete[] rgb;
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
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	double box[5][5] = { {1, 2, 3, 2, 1},
						 {2, 4, 6, 4, 2},
						 {3, 6, 9, 6, 3},
						 {2, 4, 6, 4, 2},
						 {1, 2, 3, 2, 1} };
	for (int j = 0; j < 5; j++)
		for (int i = 0; i < 5; i++)
			box[j][i] = (box[j][i] / 81.0);
	unsigned char* rgb = To_RGB();
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width; y++) {
			for (int color = 0; color < 3; color++) {
				double total = 0;
				for (int r = 0; r < 5; r++) {
					for (int c = 0; c < 5; c++) {
						int rPos = x - 2 + r;
						int cPos = y - 2 + c;
						if (rPos < 0) {
							rPos = -rPos;
						}
						if (rPos >= height) {
							rPos = ((height * 2) - 2) - rPos;
						}
						if (cPos < 0) {
							cPos = -cPos;
						}
						if (cPos >= width) {
							cPos = ((width * 2) - 2) - cPos;
						}
						total += rgb[((rPos * width + cPos) * 3) + color] * box[r][c];
					}
				}
				data[(x * width + y) * 4 + color] = (unsigned char)total;
				if (total > 255) {
					data[(x * width + y) * 4 + color] = 255;
				}
				else if (total < 0) {
					data[(x * width + y) * 4 + color] = 0;
				}
			}
		}
	}
	delete[] rgb;
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
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	double box[5][5] = { {1.0, 4.0, 6.0, 4.0, 1.0},
					   {4.0, 16.0, 24.0, 16.0, 4.0},
					   {6.0, 24.0, 36.0, 24.0, 6.0},
					   {4.0, 16.0, 24.0, 16.0, 4.0},
					   {1.0, 4.0, 6.0, 4.0, 1.0} };
	for (int j = 0; j < 5; j++)
		for (int i = 0; i < 5; i++)
			box[j][i] = (box[j][i] / 256.0);
	unsigned char* rgb = To_RGB();
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width; y++) {
			for (int color = 0; color < 3; color++) {
				double total = 0;
				for (int r = 0; r < 5; r++) {
					for (int c = 0; c < 5; c++) {
						int rPos = x - 2 + r;
						int cPos = y - 2 + c;
						if (rPos < 0) {
							rPos = -rPos;
						}
						if (rPos >= height) {
							rPos = ((height * 2) - 2) - rPos;
						}
						if (cPos < 0) {
							cPos = -cPos;
						}
						if (cPos >= width) {
							cPos = ((width * 2) - 2) - cPos;
						}
						total += rgb[((rPos * width + cPos) * 3) + color] * box[r][c];
					}
				}
				data[(x * width + y) * 4 + color] = (unsigned char)total;
				if (total > 255) {
					data[(x * width + y) * 4 + color] = 255;
				}
				else if (total < 0) {
					data[(x * width + y) * 4 + color] = 0;
				}
			}
		}
	}
	delete[] rgb;
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
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	std::vector<double> gauss;
	gauss.resize(N);
	for (size_t n = 0; n < N; n++) {
		gauss[n] = (double)(Binomial(N - 1, n));
	}

	std::vector< std::vector<double> > box;//Easy dynamic 2D filter
	box.resize(N);
	for (size_t i = 0; i < N; i++) {
		box[i].resize(N);
	}

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			box[i][j] = gauss[i] * gauss[j];
		}
	}
	for (size_t j = 0; j < N; j++)
		for (size_t i = 0; i < N; i++)
			box[j][i] = (box[j][i] / pow(16, (N / 2)));

	unsigned char* rgb = To_RGB();
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width; y++) {
			for (int color = 0; color < 3; color++) {
				double total = 0;
				for (size_t r = 0; r < N; r++) {
					for (size_t c = 0; c < N; c++) {
						int rPos = x - (N / 2) + r;
						int cPos = y - (N / 2) + c;
						if (rPos < 0) {
							rPos = -rPos;
						}
						if (rPos >= height) {
							rPos = ((height * 2) - 2) - rPos;
						}
						if (cPos < 0) {
							cPos = -cPos;
						}
						if (cPos >= width) {
							cPos = ((width * 2) - 2) - cPos;
						}
						total += rgb[((rPos * width + cPos) * 3) + color] * box[r][c];
					}
				}
				data[(x * width + y) * 4 + color] = (unsigned char)total;
				if (total > 255) {
					data[(x * width + y) * 4 + color] = 255;
				}
				else if (total < 0) {
					data[(x * width + y) * 4 + color] = 0;
				}
			}
		}
	}
	delete[] rgb;
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
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	double box[5][5] = { {-1.0, -4.0, -6.0, -4.0, -1.0},
					   {-4.0, -16.0, -24.0, -16.0, -4.0},
					   {-6.0, -24.0, 220.0, -24.0, -6.0},
					   {-4.0, -16.0, -24.0, -16.0, -4.0},
					   {-1.0, -4.0, -6.0, -4.0, -1.0} };
	for (int j = 0; j < 5; j++)
		for (int i = 0; i < 5; i++)
			box[j][i] = (box[j][i] / 256.0);
	unsigned char* rgb = To_RGB();
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width; y++) {
			for (int color = 0; color < 3; color++) {
				double total = 0;
				for (int r = 0; r < 5; r++) {
					for (int c = 0; c < 5; c++) {
						int rPos = x - 2 + r;
						int cPos = y - 2 + c;
						if (rPos < 0) {
							rPos = -rPos;
						}
						if (rPos >= height) {
							rPos = ((height * 2) - 2) - rPos;
						}
						if (cPos < 0) {
							cPos = -cPos;
						}
						if (cPos >= width) {
							cPos = ((width * 2) - 2) - cPos;
						}
						total += rgb[((rPos * width + cPos) * 3) + color] * box[r][c];
					}
				}
				data[(x * width + y) * 4 + color] = (unsigned char)total;
				if (total > 255) {
					data[(x * width + y) * 4 + color] = 255;
				}
				else if (total < 0) {
					data[(x * width + y) * 4 + color] = 0;
				}
			}
		}
	}
	delete[] rgb;
	return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	double box[5][5] = { {-1.0, -4.0, -6.0, -4.0, -1.0},
				   {-4.0, -16.0, -24.0, -16.0, -4.0},
				   {-6.0, -24.0, 476.0, -24.0, -6.0},
				   {-4.0, -16.0, -24.0, -16.0, -4.0},
				   {-1.0, -4.0, -6.0, -4.0, -1.0} };
	for (int j = 0; j < 5; j++)
		for (int i = 0; i < 5; i++)
			box[j][i] = (box[j][i] / 256.0);
	unsigned char* rgb = To_RGB();
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width; y++) {
			for (int color = 0; color < 3; color++) {
				double total = 0;
				for (int r = 0; r < 5; r++) {
					for (int c = 0; c < 5; c++) {
						int rPos = x - 2 + r;
						int cPos = y - 2 + c;
						if (rPos < 0) {
							rPos = -rPos;
						}
						if (rPos >= height) {
							rPos = ((height * 2) - 2) - rPos;
						}
						if (cPos < 0) {
							cPos = -cPos;
						}
						if (cPos >= width) {
							cPos = ((width * 2) - 2) - cPos;
						}
						total += rgb[((rPos * width + cPos) * 3) + color] * box[r][c];
					}
				}
				data[(x * width + y) * 4 + color] = (unsigned char)total;
				if (total > 255) {
					data[(x * width + y) * 4 + color] = 255;
				}
				else if (total < 0) {
					data[(x * width + y) * 4 + color] = 0;
				}
			}
		}
	}
	delete[] rgb;
	return true;
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
	if (!data)
	{
		cout << "No data." << endl;
		return NULL;
	}
	unsigned char* copy = new unsigned char[width * height];
	double box[3][3] = { {1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0},
						{1.0 / 8.0, 1.0 / 4.0, 1.0 / 8.0},
						{1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0} };
	unsigned char* rgb = To_RGB();
	for (int x = 0; x < height; x += 2) {
		for (int y = 0; y < width; y += 2) {
			for (int color = 0; color < 3; color++) {
				double total = 0;
				for (int r = 0; r < 3; r++) {
					for (int c = 0; c < 3; c++) {
						int rPos = x - 1 + r;
						int cPos = y - 1 + c;
						if (rPos < 0) {
							rPos = -rPos;
						}
						if (rPos >= height) {
							rPos = ((height * 2) - 1) - rPos;
						}
						if (cPos < 0) {
							cPos = -cPos;
						}
						if (cPos >= width) {
							cPos = ((width * 2) - 1) - cPos;
						}
						total += rgb[((rPos * width + cPos) * 3) + color] * box[r][c];
					}
				}
				copy[(x / 2 * width / 2 + y / 2) * 4 + color] = (unsigned char)total;
			}
			copy[(x / 2 * width / 2 + y / 2) * 4 + 4] = 255;
		}
	}
	height = height / 2;
	width = width / 2;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			unsigned char* dest = data + ((y * width + x) * 4);
			for (int i = 0; i < 3; i++) {
				dest[i] = copy[(y * width + x) * 4 + i];
			}

		}
	}
	delete[] rgb;
	return true;

}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
	ClearToBlack();
	return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
	ClearToBlack();
	return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
	ClearToBlack();
	return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
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

