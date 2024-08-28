#ifndef JPEGHANDLER_H
#define JPEGHANDLER_H
#include <stdio.h>
#include <jpeglib.h>
#ifndef TRUE
#define TRUE 1
#endif // !TRUE
#ifndef FALSE
#define FALSE 0
#endif // !FALSE
/*returns the image dimensions and color space*/
int
jpegh_isjpeg(FILE *fp);
int
jpegh_get_image_info(FILE* in_file, size_t* image_width, size_t* image_height,  
	int* num_components, int* color_space);

/*writes image data to a jpeg file*/
int
jpegh_write_image(FILE* out_file, unsigned char** img_data,  size_t image_width, 
	size_t image_height, int num_components);

/*gets image data from a jpeg data stream*/
unsigned char**
jpegh_get_image(FILE* in_file, size_t* image_width, size_t* image_height, 
	int* num_components);

#endif // !READIMGINFO_H