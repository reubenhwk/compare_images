
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <magick/api.h>

#include "fileio.h"

/**
 * @brief read an image file into a 32bit float array, converted to gray
 */
float *read_img_f32_gray(const char *fname, size_t * nx, size_t * ny)
{

	ExceptionInfo *exception = AcquireExceptionInfo();
	ImageInfo *image_info = CloneImageInfo((ImageInfo *) NULL);
	unsigned char *file_data = 0;
	size_t file_data_size = read_file(&file_data, fname);
	Image *image =
	    BlobToImage(image_info, file_data, file_data_size, exception);

	if (exception->severity != UndefinedException)
		CatchException(exception);
	if (!image) {
		fprintf(stderr, "failed to read image to a blob.\n");
		exit(-1);
	}

	if (file_data)
		free(file_data);

	unsigned char *rgb =
	    malloc(3 * image->magick_columns * image->magick_rows);
	ExportImagePixels(image, 0, 0, image->magick_columns,
			  image->magick_rows, "RGB", CharPixel, rgb, exception);

	float *retval =
	    malloc(sizeof(float) * image->magick_columns * image->magick_rows);
#pragma omp parallel for
	for (size_t i = 0; i < image->magick_columns * image->magick_rows; ++i) {
		retval[i] = (6969.0 * rgb[3 * i + 0]
			     + 23434.0 * rgb[3 * i + 1]
			     + 2365.0 * rgb[3 * i + 2]) / 32768.0;
	}

	if (exception->severity != UndefinedException)
		CatchException(exception);

	*nx = image->magick_columns;
	*ny = image->magick_rows;

	image = DestroyImage(image);
	image_info = DestroyImageInfo(image_info);
	exception = DestroyExceptionInfo(exception);

	return retval;
}
