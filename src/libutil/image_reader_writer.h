#ifndef IMAGE_READER_WRITER_H
#define IMAGE_READER_WRITER_H
#define IMG_TYPE_PNG 0
#define IMG_TYPE_JPEG 1
unsigned char **
img_get_image(char *file_name, size_t *width, size_t *height,
int *num_of_channels);
int 
img_write_image(unsigned char **image,char *file_name, size_t width, 
size_t height, int num_of_channels, int type);
void
img_free_image(void **image);
double**
img_char_to_double(unsigned char **ch_image, size_t width, size_t height,
int num_of_channels);
unsigned char**
img_double_to_char(double** do_image, size_t width, size_t height,
int num_of_channels);
void **
img_get_row_pointers(void *compact_img,int width, int height,int nchannels, int nbytes);
#endif