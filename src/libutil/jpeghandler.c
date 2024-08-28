#include <stdlib.h>
#include "jpeghandler.h"
#include "image_reader_writer.h"
#define JPEGH_MAX_IMAGE_HEIGHT 5000
#define JPEGH_MAX_IMAGE_WIDTH 5000
/*checks whether a give file is a jpeg data stream*/
int
jpegh_isjpeg(FILE *fp)
{
	int is_jpeg = 0, i = 0;
	size_t file_size;
	unsigned char start_maker[2];
	if(fp == NULL)
	{
		return (0);
	}
	fseek(fp,0L,SEEK_END);
	file_size = ftell(fp);
	fseek(fp,0L, SEEK_SET);
	/*search for start of image imaker*/
	while(i < file_size)
	{
		fread(start_maker, sizeof(char),2,fp);
		is_jpeg = (start_maker[0] == 0xFF) && (start_maker[1] == 0xD8);
		if(is_jpeg)
		{
			break;
		}
		i += 2;
	}
	fseek(fp,0L, SEEK_SET);
	return (is_jpeg);
}

/*returns the image dimensions and color space*/
int
jpegh_get_image_info(FILE* in_file, size_t* image_width, size_t* image_height,  
	int* num_components, int* color_space)
{
	if (in_file == NULL)
	{
		return (FALSE); /*operation failed*/
	}
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, in_file);
	jpeg_read_header(&cinfo, TRUE);
	*image_height = cinfo.image_height;
	*image_width = cinfo.image_width;
	*num_components = cinfo.num_components;
	*color_space = cinfo.jpeg_color_space;
	jpeg_destroy_decompress(&cinfo);
	return (TRUE);
}

/*writes image data to a jpeg file*/
int
jpegh_write_image(FILE* out_file, unsigned char** img_data,  size_t image_width, 
	size_t image_height, int num_components)
{
	if (out_file == NULL)
	{
		return (FALSE); /*process failed*/
	}
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, out_file);
	cinfo.image_height = image_height;
	cinfo.image_width = image_width;
	cinfo.input_components = num_components;
	switch (num_components)
	{
	case 1:
	cinfo.in_color_space = JCS_GRAYSCALE;
		break;
	case 3:
		cinfo.in_color_space = JCS_RGB;
		break;
	default:
		break;
	}
	jpeg_set_defaults(&cinfo);
	jpeg_start_compress(&cinfo, TRUE);
	JSAMPROW row_pointer[1];
	size_t i = 0;
	while (cinfo.next_scanline < cinfo.image_height && i < cinfo.image_height) 
	{
		row_pointer[0] = (JSAMPROW)img_data[i++];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	return (TRUE);
}

/*gets image data from a jpeg data stream*/
unsigned char**
jpegh_get_image(FILE* in_file, size_t* image_width, size_t* image_height, 
	int* num_components) 
	{
	if (in_file == NULL)
	{
		return (NULL);
	}
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, in_file);
	jpeg_read_header(&cinfo, TRUE);
	J_COLOR_SPACE color_space;
	jpeg_start_decompress(&cinfo);
	*image_height = cinfo.output_height;
	*image_width = cinfo.output_width;
	*num_components = cinfo.output_components;
	color_space = cinfo.out_color_space;
	int exit = (cinfo.output_width > JPEGH_MAX_IMAGE_WIDTH) || 
	(cinfo.output_height > JPEGH_MAX_IMAGE_HEIGHT) || 
	(color_space != JCS_RGB && color_space != JCS_GRAYSCALE);
	if (exit)
	{
		jpeg_finish_decompress(&cinfo);
		jpeg_destroy_decompress(&cinfo);
		return (NULL);
	}
	size_t row_stride = cinfo.output_width * cinfo.output_components * sizeof(char);

	unsigned char** img_data = 
	(unsigned char**)img_get_row_pointers(malloc(row_stride * cinfo.output_height),
											cinfo.output_width,cinfo.output_height,
											cinfo.output_components,sizeof(char));
	if(!img_data)
	{
		jpeg_finish_decompress(&cinfo);
		jpeg_destroy_decompress(&cinfo);
		return (NULL);
	}
	// TODO: add functionality to read exif data, and rotate image accordinly
	JSAMPROW row_pointer[1];
	size_t i = 0;
	while (cinfo.output_scanline < cinfo.output_height && i < cinfo.output_height)
	{
		row_pointer[0] = (JSAMPROW)img_data[i++];
		jpeg_read_scanlines(&cinfo,row_pointer, 1);
	}
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	return (img_data);
}
