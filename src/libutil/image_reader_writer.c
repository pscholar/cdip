#include <stdio.h>
#include <stdlib.h>
#include "image_reader_writer.h"
#include "jpeghandler.h"
#include "pnghandler.h"

unsigned char **
img_get_image(char *file_name, size_t *width, size_t *height,
int *num_of_channels)
{
  unsigned char **image = NULL;
  FILE *fp;
  fp = fopen(file_name,"rb");
  if(fp == NULL)
  {
    return NULL;
  }
  if(pngh_ispng(fp))
  {
    image = pngh_get_image(fp,width,height,num_of_channels);
  }
  else if (jpegh_isjpeg(fp))
  {
    image = jpegh_get_image(fp,width,height,num_of_channels);
  }
  else
  {
    printf("unsupported file format\n");
  } 
  fclose(fp);
  return (image);
}

int 
img_write_image(unsigned char **image,char *file_name, size_t width, 
size_t height,int num_of_channels, int type)
{
  FILE *fp = fopen(file_name,"wb");
  int  success = 0;
  if(fp == NULL)
  {
    return (success);
  }
  if(type == IMG_TYPE_PNG)
  {
    pngh_write_image(fp,image,width,height, num_of_channels);
  }else if(type == IMG_TYPE_JPEG)
  {
    jpegh_write_image(fp,image,width,height,num_of_channels);
  }
  fclose(fp);
  return (success);
}


void
img_free_image(void **image)
{
  if(image == NULL) {return;};
  free(image[0]);
  free(image);
}


double**
img_char_to_double(unsigned char **ch_image, size_t width, size_t height,
int num_of_channels)
{
  size_t row_bytes = ((width) * (num_of_channels));
  double **do_image = (double **)malloc(sizeof(double *) * (height));
  for (size_t i = 0; i < height; i++)
  {
    do_image[i] = (double *)malloc(sizeof(double) * row_bytes);
    for(size_t j = 0; j < row_bytes ; j++)
    {
      do_image[i][j] = (double)ch_image[i][j];
    }
  }
  return (do_image);
}


unsigned char**
img_double_to_char(double **do_image, size_t width, size_t height,
int num_of_channels)
{
  size_t row_bytes = ((width)  * (num_of_channels));
  unsigned char **ch_image = (unsigned char **)malloc(sizeof(char *) * (height));
  for (size_t i = 0; i < height; i++)
  {
    ch_image[i] = (unsigned char *)malloc(sizeof(char) * row_bytes);
    for(size_t j = 0; j < row_bytes ; j++)
    {
      ch_image[i][j] = (unsigned char)do_image[i][j];
    }
  }
  return (ch_image);
}

void **
img_get_row_pointers(void *compact_img, int width, int height, int nchannels, int nbytes)
{
  if(!compact_img){ return (NULL);}

  unsigned char **img = malloc(sizeof(void *) * height);
  if(!img){ return (NULL);}

  size_t warp;
  warp = width * nchannels * nbytes;

  for (size_t i = 0; i < height; i++)
  {
    img[i] = (unsigned char *)(compact_img + i * warp);
  }

  return (void **)img;
}