
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

float *read_img_f32_gray(const char *fname, size_t * nx, size_t * ny);
unsigned char *read_img_rgb(const char *fname, size_t * nx, size_t * ny);
unsigned char *read_img_rgba(const char *fname, size_t * nx, size_t * ny);

#ifdef __cplusplus
}
#endif
