// Copyright (c) 2008-2011, Guoshen Yu <yu@cmap.polytechnique.fr>
// Copyright (c) 2008-2011, Jean-Michel Morel <morel@cmla.ens-cachan.fr>
//
// WARNING:
// This file implements an algorithm possibly linked to the patent
//
// Jean-Michel Morel and Guoshen Yu, Method and device for the invariant
// affine recognition recognition of shapes (WO/2009/150361), patent pending.
//
// This file is made available for the exclusive aim of serving as
// scientific tool to verify of the soundness and
// completeness of the algorithm description. Compilation,
// execution and redistribution of this file may violate exclusive
// patents rights in certain countries.
// The situation being different for every country and changing
// over time, it is your responsibility to determine which patent
// rights restrictions apply to you before you compile, use,
// modify, or redistribute this file. A patent lawyer is qualified
// to make this determination.
// If and only if they don't conflict with any patent terms, you
// can benefit from the following license terms attached to this
// file.
//
// This program is provided for scientific and educational only:
// you can use and/or modify it for these purposes, but you are
// not allowed to redistribute this work or derivative works in
// source or executable form. A license must be obtained from the
// patent right holders for any other use.
//
//
//*------------------------ compute_asift_keypoints -------------------------*/
// Compute the ASIFT keypoints on the input image.
//
// Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
//
// Reference: J.M. Morel and G.Yu, ASIFT: A New Framework for Fully Affine Invariant Image
//            Comparison, SIAM Journal on Imaging Sciences, vol. 2, issue 2, pp. 438-469, 2009.
// Reference: ASIFT online demo (You can try ASIFT with your own images online.)
//                        http://www.ipol.im/pub/algo/my_affine_sift/
/*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "compute_asift_keypoints.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define ABS(x)    (((x) > 0) ? (x) : (-(x)))

/* InitSigma gives the amount of smoothing applied to the image at the
first level of each octave.  In effect, this determines the sampling
needed in the image domain relative to amount of smoothing.  Good
values determined experimentally are in the range 1.2 to 1.8.
*/
/* float InitSigma_aa = 1.0;*/
static float InitSigma_aa = 1.6;

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

/* Gaussian convolution kernels are truncated at this many sigmas from
the center.  While it is more efficient to keep this value small,
experiments show that for consistent scale-space analysis it needs
a value of about 3.0, at which point the Gaussian has fallen to
only 1% of its central value.  A value of 2.0 greatly reduces
keypoint consistency, and a value of 4.0 is better than 3.0.
*/
static const float GaussTruncate1 = 4.0;

/* --------------------------- Blur image --------------------------- */

/* Same as ConvBuffer, but implemented with loop unrolling for increased
speed.  This is the most time intensive routine in keypoint detection,
so deserves careful attention to efficiency.  Loop unrolling simply
sums 5 multiplications at a time to allow the compiler to schedule
operations better and avoid loop overhead.  This almost triples
speed of previous version on a Pentium with gcc.
*/
static inline void ConvBufferFast(float *buffer, float const * const kernel, int rsize, int ksize)
{
	for (int i = 0; i < rsize; ++i) {
		buffer[i] = 0.0;
		for (int j = 0; j < ksize; ++j) {
			buffer[i] += buffer[i+j] * kernel[j];
		}
	}
}

/* Convolve image with the 1-D kernel vector along image rows.  This
is designed to be as efficient as possible.  Pixels outside the
image are set to the value of the closest image pixel.
*/
static void ConvHorizontal(vector < float >&image, int width, int height, float const * const kernel, int ksize)
{
	int const rows = height;
	int const cols = width;

	int const halfsize = ksize / 2;

#pragma omp parallel for
	for (int r = 0; r < rows; r++) {
		vector<float>buffer;
		buffer.resize(ksize + cols);

		/* Copy the row into buffer with pixels at ends replicated for
		   half the mask size.  This avoids need to check for ends
		   within inner loop. */
		for (int i = 0; i < halfsize; i++) {
			buffer[i] = image[r * cols];
			buffer[halfsize + cols + i] = image[r * cols + cols - 1];
		}
		memcpy(&buffer[halfsize], &image[r * cols], sizeof(float) * cols);
		ConvBufferFast(&buffer[0], kernel, cols, ksize);
		memcpy(&image[r * cols], &buffer[0], sizeof(float) * cols);
	}
}

/* Same as ConvHorizontal, but apply to vertical columns of image.
*/
static void ConvVertical(vector < float >&image, int width, int height, float const * const kernel, int ksize)
{
	int const rows = height;
	int const cols = width;
	int const halfsize = ksize / 2;

#pragma omp parallel for
	for (int c = 0; c < cols; c++) {
		vector<float>buffer;
		buffer.resize(ksize + rows);

		for (int i = 0; i < halfsize; i++) {
			buffer[i] = image[c];
			buffer[halfsize + rows + i] = image[(rows - 1) * cols + c];
		}
		for (int i = 0; i < rows; i++)
			buffer[halfsize + i] = image[i * cols + c];

		ConvBufferFast(&buffer[0], kernel, rows, ksize);

		for (int r = 0; r < rows; r++) {
			image[r * cols + c] = buffer[r];
		}
	}
}

/* 1D Convolve image with a Gaussian of width sigma and store result back
in image.   This routine creates the Gaussian kernel, and then applies
it in horizontal (flag_dir=0) OR vertical directions (flag_dir!=0).
*/
void GaussianBlur1D(vector < float >&image, int width, int height, float sigma, int flag_dir)
{
	float x, kernel[100], sum = 0.0;
	int ksize, i;

	/* The Gaussian kernel is truncated at GaussTruncate sigmas from
	   center.  The kernel size should be odd.
	 */
	ksize = (int)(2.0 * GaussTruncate1 * sigma + 1.0);
	ksize = MAX(3, ksize);	/* Kernel must be at least 3. */
	ksize = MIN(99, ksize);	/* Kernel must be no more than 100. */
	if (ksize % 2 == 0)	/* Make kernel size odd. */
		ksize++;
	assert(ksize < 100);

	/* Fill in kernel values. */
	for (i = 0; i <= ksize; i++) {
		x = i - ksize / 2;
		kernel[i] = exp(-x * x / (2.0 * sigma * sigma));
		sum += kernel[i];
	}
	/* Normalize kernel values to sum to 1.0. */
	for (i = 0; i < ksize; i++)
		kernel[i] /= sum;

	if (flag_dir == 0) {
		ConvHorizontal(image, width, height, kernel, ksize);
	} else {
		ConvVertical(image, width, height, kernel, ksize);
	}
}

static inline void compensate_affine_coor1(float *x0, float *y0, int w1, int h1, float t1, float t2, float Rtheta)
{
	float x_ori, y_ori;

	float x1 = *x0;
	float y1 = *y0;

	Rtheta = Rtheta * PI / 180.0f;
	float const sin_Rtheta = sin(Rtheta);
	float const cos_Rtheta = cos(Rtheta);

	if (Rtheta <= PI / 2.0f) {
		x_ori = 0.0f;
		y_ori = w1 * sin_Rtheta / t1;
	} else {
		x_ori = -w1 * cos_Rtheta / t2;
		y_ori = (w1 * sin_Rtheta + h1 * sin(Rtheta - PI / 2.0f)) / t1;
	}

	/* project the coordinates of im1 to original image before tilt-rotation transform */
	/* Get the coordinates with respect to the 'origin' of the original image before transform */
	x1 = x1 - x_ori;
	y1 = y1 - y_ori;
	/* Invert tilt */
	x1 = x1 * t2;
	y1 = y1 * t1;
	/* Invert rotation (Note that the y direction (vertical) is inverse to the usual concention. Hence Rtheta instead of -Rtheta to inverse the rotation.) */
	*x0 = cos_Rtheta * x1 - sin_Rtheta * y1;
	*y0 = sin_Rtheta * x1 + cos_Rtheta * y1;
}

/* -------------- MAIN FUNCTION ---------------------- */

KeyPoints compute_asift_keypoints(vector < float > const &image, int width, int height, int verb,
			    siftPar & siftparameters)
// Compute ASIFT keypoints in the input image.
// Input:
// image: input image
// width, height: width and height of the input image.
// num_of_tilts: number of tilts to simulate.
// verb: 1/0 --> show/don not show verbose messages. (1 for debugging)
// keys_all (output): ASIFT keypoints. It is a 2D matrix with varying rows and columns. Each entry keys_all[tt][rr]
//      stores the SIFT keypoints calculated on the image with the simulated tilt index tt and simulated rotation index rr (see the code below). In the coordinates of the keypoints,
//      the affine distortions have been compensated.
// siftparameters: SIFT parameters.
//
// Output: the number of keypoints
{
	// number N of tilts to simulate t = 1, \sqrt{2}, (\sqrt{2})^2, ..., {\sqrt{2}}^(N-1)
	int const num_of_tilts = 7;
	vector < vector < keypointslist > >keys_all;
	float t_min, t_k;
	int num_tilt, num_rot_t2;
	int fproj_o;
	float fproj_p, fproj_bg;
	char fproj_i;
	float *fproj_x4, *fproj_y4;
	//  float frot_b=0;
	float frot_b = 128;
	char *frot_k;
	int counter_sim = 0, num_sim;
	int flag_dir = 1;
	float BorderFact = 6 * sqrt(2.);


	fproj_o = 3;
	fproj_p = 0;
	fproj_i = 0;
	fproj_bg = 0;
	fproj_x4 = 0;
	fproj_y4 = 0;

	frot_k = 0;

	num_rot_t2 = 10;

	t_min = 1;
	t_k = sqrt(2.);

	num_tilt = num_of_tilts;

	if (num_tilt < 1) {
		printf("Number of tilts num_tilt should be equal or larger than 1. \n");
		exit(-1);
	}

	/* Calculate the number of simulations, and initialize keys_all */
	keys_all = std::vector < vector < keypointslist > >(num_tilt);
#pragma omp parallel for
	for (int tt = 1; tt <= num_tilt; tt++) {
		float t = t_min * pow(t_k, tt - 1);

		if (t == 1) {
			counter_sim++;

			keys_all[tt - 1] = std::vector < keypointslist > (1);
		} else {
			int num_rot1 = round(num_rot_t2 * t / 2);
			if (num_rot1 % 2 == 1) {
				num_rot1 = num_rot1 + 1;
			}
			num_rot1 = num_rot1 / 2;
			counter_sim += num_rot1;

			keys_all[tt - 1] = std::vector < keypointslist > (num_rot1);
		}
	}

	num_sim = counter_sim;

	if (verb) {
		printf("%d affine simulations will be performed. \n", num_sim);
	}

	counter_sim = 0;

	/* Affine simulation (rotation+tilt simulation) */
	// Loop on tilts.
#ifdef _OPENMP
	omp_set_nested(1);
#endif
#pragma omp parallel for
	for (int tt = 1; tt <= num_tilt; tt++) {
		float t = t_min * pow(t_k, tt - 1);

		float t1 = 1;
		float t2 = 1 / t;

		// If tilt t = 1, do not simulate rotation.
		if (t == 1) {
			compute_sift_keypoints(&image[0], keys_all[tt - 1][0], width, height, siftparameters);

		} else {
			// The number of rotations to simulate under the current tilt.
			int num_rot1 = round(num_rot_t2 * t / 2);

			if (num_rot1 % 2 == 1) {
				num_rot1 = num_rot1 + 1;
			}
			num_rot1 = num_rot1 / 2;
			float delta_theta = PI / num_rot1;

			// Loop on rotations.
#pragma omp parallel for
			for (int rr = 1; rr <= num_rot1; rr++) {
				float theta = delta_theta * (rr - 1);
				theta = theta * 180 / PI;

				vector < float >image_rot;
				int width_r, height_r;

				// simulate a rotation: rotate the image with an angle theta. (the outside of the rotated image are padded with the value frot_b)
				frot(image, image_rot, width, height, &width_r, &height_r, &theta, &frot_b, frot_k);

				/* Tilt */
				int width_t = (int)(width_r * t1);
				int height_t = (int)(height_r * t2);

				int fproj_sx = width_t;
				int fproj_sy = height_t;

				float fproj_x1 = 0;
				float fproj_y1 = 0;
				float fproj_x2 = width_t;
				float fproj_y2 = 0;
				float fproj_x3 = 0;
				float fproj_y3 = height_t;

				/* Anti-aliasing filtering along vertical direction */
				/* sigma_aa = InitSigma_aa * log2(t); */
				float sigma_aa = InitSigma_aa * t / 2;
				GaussianBlur1D(image_rot, width_r, height_r, sigma_aa, flag_dir);

				// simulate a tilt: subsample the image along the vertical axis by a factor of t.
				vector < float >image_tilt(width_t * height_t);
				fproj(image_rot, image_tilt, width_r, height_r, &fproj_sx, &fproj_sy, &fproj_bg, &fproj_o,
				      &fproj_p, &fproj_i, fproj_x1, fproj_y1, fproj_x2, fproj_y2, fproj_x3, fproj_y3,
				      fproj_x4, fproj_y4);

				if (verb) {
					printf("Rotation theta = %.2f, Tilt t = %.2f. w=%d, h=%d, sigma_aa=%.2f, \n", theta,
					       t, width_t, height_t, sigma_aa);
				}

				// compute SIFT keypoints_unfiltered on simulated image.
				keypointslist keypoints_unfiltered;
				compute_sift_keypoints(&image_tilt[0], keypoints_unfiltered, width_t, height_t, siftparameters);

				/* check if the keypoint is located on the boundary of the parallelogram (i.e., the boundary of the distorted input image). If so, remove it to avoid boundary artifacts. */
				if (keypoints_unfiltered.size() != 0) {
					keypointslist keypoints_filtered;
					for (int cc = 0; cc < (int)keypoints_unfiltered.size(); cc++) {

						float x[4];
						float y[4];

						float const theta1 = theta * PI / 180;
						float const sin_theta1 = sin(theta1);
						float const cos_theta1 = cos(theta1);

						/* the coordinates of the 4 submits of the parallelogram */
						if (theta <= 90) {
							x[0] = height * sin_theta1;
							y[0] = 0;
							y[1] = width * sin_theta1;
							x[2] = width * cos_theta1;
							x[3] = 0;
							y[3] = height * cos_theta1;
							x[1] = x[0] + x[2];
							y[2] = y[1] + y[3];

							/* note that the vertical direction goes from top to bottom!!!
							   The calculation above assumes that the vertical direction goes from the bottom to top. Thus the vertical coordinates need to be reversed!!! */
							y[0] = y[2] - y[0];
							y[1] = y[2] - y[1];
							y[3] = y[2] - y[3];
							y[2] = 0;

							y[0] = y[0] * t2;
							y[1] = y[1] * t2;
							y[2] = y[2] * t2;
							y[3] = y[3] * t2;
						} else {
							y[0] = -height * cos_theta1;
							x[1] = height * sin_theta1;
							x[2] = 0;
							y[2] = width * sin_theta1;
							x[3] = -width * cos_theta1;
							y[3] = 0;
							x[0] = x[1] + x[3];
							y[1] = y[0] + y[2];

							/* note that the vertical direction goes from top to bottom!!!
							   The calculation above assumes that the vertical direction goes from the bottom to top. Thus the vertical coordinates need to be reversed!!! */
							y[0] = y[1] - y[0];
							y[2] = y[1] - y[2];
							y[3] = y[1] - y[3];
							y[1] = 0;

							y[0] = y[0] * t2;
							y[1] = y[1] * t2;
							y[2] = y[2] * t2;
							y[3] = y[3] * t2;
						}

						float x0 = keypoints_unfiltered[cc].x;
						float y0 = keypoints_unfiltered[cc].y;

						/* the distances from the keypoint to the 4 sides of the parallelogram */
						float const d1 = ABS((x[1] - x[0]) * (y[0] - y0) -
							 (x[0] - x0) * (y[1] - y[0])) / sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] -
														y[0]) * (y[1] -
														       y[0]));
						float const d2 = ABS((x[2] - x[1]) * (y[1] - y0) -
							 (x[1] - x0) * (y[2] - y[1])) / sqrt((x[2] - x[1]) * (x[2] - x[1]) + (y[2] -
														y[1]) * (y[2] -
														       y[1]));
						float const d3 = ABS((x[3] - x[2]) * (y[2] - y0) -
							 (x[2] - x0) * (y[3] - y[2])) / sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] -
														y[2]) * (y[3] -
														       y[2]));
						float const d4 = ABS((x[0] - x[3]) * (y[3] - y0) -
							 (x[3] - x0) * (y[0] - y[3])) / sqrt((x[0] - x[3]) * (x[0] - x[3]) + (y[0] -
														y[3]) * (y[0] -
														       y[3]));

						float const scale1 = keypoints_unfiltered[cc].scale;
						float const BorderTh = BorderFact * scale1;

						if (!
						    ((d1 < BorderTh) || (d2 < BorderTh) || (d3 < BorderTh)
						     || (d4 < BorderTh))) {
							// Normalize the coordinates of the matched points by compensate the simulate affine transformations
							compensate_affine_coor1(&x0, &y0, width, height, 1 / t2, t1, theta);
							keypoints_unfiltered[cc].x = x0;
							keypoints_unfiltered[cc].y = y0;

							keypoints_filtered.push_back(keypoints_unfiltered[cc]);
						}
					}
					keys_all[tt - 1][rr - 1] = keypoints_filtered;
				}
			}
		}
	}


	int num_keys_total = 0;
#pragma omp parallel for reduction(+:num_keys_total)
	for (int tt = 0; tt < (int)keys_all.size(); tt++) {
		for (int rr = 0; rr < (int)keys_all[tt].size(); rr++) {
			num_keys_total += (int)keys_all[tt][rr].size();
		}
	}
	printf("%d ASIFT keypoints are detected. \n", num_keys_total);

	KeyPoints retval;
	retval.keys = keys_all;
	retval.tilts = num_of_tilts;
	retval.count = num_keys_total;
	return retval;
}
