#include <stdio.h>
#include <stdlib.h>
#include <png.h>
#include "pnghandler.h"
#define IMAGE_MAX_SUPPORTED_HEIGHT 5000
#define IMAGE_MAX_SUPPORTED_WIDTH 5000
#define NUM_PNG_SIGNATURE_BYTES 8

/*checks whether a given file is a png datastream file*/
int
pngh_ispng(FILE *fp)
{
  unsigned char header[NUM_PNG_SIGNATURE_BYTES];
  int is_png;
  if(fp == NULL)
  {
    return (0);
  }
  /*move back to the start of the file*/
  fseek(fp,0L, SEEK_SET);
  if( fread(header,sizeof(char),NUM_PNG_SIGNATURE_BYTES, fp) 
  != NUM_PNG_SIGNATURE_BYTES )
  {
    return (0);
  }
  is_png = (png_sig_cmp(header, 0, NUM_PNG_SIGNATURE_BYTES) == 0);
  /*move back to the start of the file*/
  fseek(fp,0L, SEEK_SET);
  return (is_png);
}

int
pngh_get_image_info(FILE *fp, size_t *width, size_t *height,
int *num_of_channels)
{
  int is_png = pngh_ispng(fp);
  if(!is_png)
  {
    return (0);
  }
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL,
  NULL, NULL);
  if (!png_ptr)
  {
    return (0);
  }
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  {
    png_destroy_read_struct(&png_ptr,(png_infopp)NULL, (png_infopp)NULL);
    //error ocurred
    return (0);
  }
  png_init_io(png_ptr, fp);
  png_set_user_limits(png_ptr, IMAGE_MAX_SUPPORTED_WIDTH, 
  IMAGE_MAX_SUPPORTED_HEIGHT); 
  png_read_info(png_ptr, info_ptr); 
  int  bit_depth, color_type, interlace_type,compression_type, 
  filter_method,channels;
  png_get_IHDR(png_ptr, info_ptr, (png_uint_32 *)width, (png_uint_32 *)height,
  &bit_depth, &color_type, &interlace_type,&compression_type, &filter_method);
  channels = png_get_channels(png_ptr, info_ptr);
  (*num_of_channels) = channels;
  png_destroy_read_struct(&png_ptr,&info_ptr, (png_infopp)NULL);
  return (1);
}

/**
 * reads image data from a png file
*/ 
unsigned char **
pngh_get_image(FILE *fp, size_t *width, size_t *height,int *num_of_channels)
{
  /*check whether the file is a png datastream file*/
  int is_png = pngh_ispng(fp);
  if(!is_png)
  {
    return (NULL);
  }
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL,
  NULL, NULL);
  if (!png_ptr)
  {
    return NULL;
  }
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  {
    png_destroy_read_struct(&png_ptr,(png_infopp)NULL, (png_infopp)NULL);
    //error ocurred
    return NULL;
  }
  png_init_io(png_ptr, fp);//set up the input code
  png_set_user_limits(png_ptr, IMAGE_MAX_SUPPORTED_WIDTH, IMAGE_MAX_SUPPORTED_HEIGHT); // set width and height limits
  png_read_info(png_ptr, info_ptr); //read file information
  int  bit_depth, color_type, interlace_type,compression_type, 
  filter_method,channels;
  png_get_IHDR(png_ptr, info_ptr, (png_uint_32 *)width, (png_uint_32 *)height,
  &bit_depth, &color_type, &interlace_type,&compression_type, &filter_method);
  channels = png_get_channels(png_ptr, info_ptr);
  if (color_type == PNG_COLOR_TYPE_PALETTE)
  {
    png_set_palette_to_rgb(png_ptr);
  }
  if (bit_depth == 16)
  {
    png_set_scale_16(png_ptr);
  }
  if (color_type & PNG_COLOR_MASK_ALPHA)
  {
    png_set_strip_alpha(png_ptr);
  }
  if (bit_depth < 8)
  {
    png_set_packing(png_ptr);
  } 
  png_set_interlace_handling(png_ptr);
  png_read_update_info(png_ptr, info_ptr);
  png_get_IHDR(png_ptr, info_ptr, (png_uint_32 *)width, (png_uint_32 *)height,
  &bit_depth, &color_type, &interlace_type,&compression_type, &filter_method);
  channels = png_get_channels(png_ptr, info_ptr);
  (*num_of_channels) = channels;
  size_t row_bytes = png_get_rowbytes(png_ptr, info_ptr); 
  unsigned char **image = (unsigned char **)malloc(sizeof(png_bytep*) * (*height));
  for (size_t i = 0; i < *height; i++)
  {
    image[i] = (unsigned char *)malloc(row_bytes);
  }
  png_read_image(png_ptr, (png_bytepp)image); 
  png_read_end(png_ptr, NULL); 
  png_destroy_read_struct(&png_ptr, &info_ptr,NULL);
  return (image); 
}

/**
 * writes png data to a png file
*/
void 
pngh_write_image(FILE *fp, unsigned char **image, size_t width, size_t height, 
int num_of_channels)
{
  if(fp == NULL)
  {
    return;
  }
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL,NULL, NULL);
  if (!png_ptr)
  {
    //error occured call error handling function
    return;
  }
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  {
    png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
    //error ocurred
    return;
  }
  png_init_io(png_ptr, fp);
  int color_type;
  if(num_of_channels == 1)
  {
    color_type = PNG_COLOR_TYPE_GRAY;
  }else
  {
    color_type = PNG_COLOR_TYPE_RGB;
  }
  png_set_IHDR(png_ptr, info_ptr, width, height,8, color_type,  PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_set_rows(png_ptr,info_ptr,(png_bytepp)image);
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fp);
  return;
}



