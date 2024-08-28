double ***image_create_image( size_t image_height, size_t image_width, int num_of_channels)
{
  double ***image = (double ***)malloc(sizeof(double **) * num_of_channels);
  for (size_t n = 0; n < num_of_channels; n++)
  {
    image[n] = image_create_image_channel(image_height,image_width);
  }
  return (image);
}
double **image_create_image_channel( size_t image_height, size_t image_width)
{
  double  **image_channel = (double **)malloc( sizeof(double *) * image_height);
  for (size_t x = 0; x < image_height; x++)
  {
    image_channel[x] = (double *)calloc( image_width,sizeof(double));
  }
  return (image_channel);
}

void image_free_image(double ***image, size_t image_height,int num_of_channels)
{
  for (size_t n = 0; n < num_of_channels; n++)
  {
    for (size_t x = 0; x < image_height; x++)
    {
      free(image[n][x]);      
    }  
  }
  free(image);
  return;
}

double ***image_rescale(double ***orig_image,size_t orig_image_height, size_t orig_image_width, size_t new_image_height, size_t new_image_width, int num_of_channels)
{
  double ***new_image = image_create_image(new_image_height,new_image_width,num_of_channels);
  int m,n;
  float horizontal_scaler = ((float)orig_image_width) / ((float)new_image_width);
  float vertical_scaler = ((float)orig_image_height) / ((float)new_image_height);
  for (int c = 0; c < num_of_channels; c++)
  {
    for (size_t x = 0; x < new_image_height; x++)
    {
      m = (int)floor(x * vertical_scaler);
      for (size_t y = 0; y < new_image_width; y++)
      {
        n = (int)floor( y * horizontal_scaler);
        new_image[c][x][y] = orig_image[c][m][n];
        
      }
    }    
  }
  return (new_image);
}