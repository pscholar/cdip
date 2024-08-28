#include <stdio.h>
#ifndef PNGHANDLER_H
#define PNGHANDLER_H
int
pngh_ispng(FILE *fp);
int
pngh_get_image_info(FILE *fp,size_t *width, size_t *height,int *num_of_channels);
unsigned char **
pngh_get_image(FILE *fp, size_t *width, size_t *height, int *num_of_channels);

void 
pngh_write_image(FILE *fp, unsigned char **image, size_t width, size_t height, 
int num_of_channels);

#endif // PNGHANDLER_H

