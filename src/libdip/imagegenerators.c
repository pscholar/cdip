double ***image_zoneplate(size_t image_height, size_t image_width, int num_of_channels, double min, double max)
{
  double ***image = image_create_image(image_height, image_width,num_of_channels);
  double range = max - min;
  double inc_x = range /(double) image_height ;
  double inc_y = range /(double) image_width ;
  double cx = image_height / 2;
  double cy = image_width / 2;
  double i,j;
  for (int n = 0; n < num_of_channels; n++)
  {
    for (size_t x = 0; x < image_height; x++)
    {
      i = ((double)x - cx) * inc_x;
      i = i * i;
      for (size_t y = 0; y < image_width; y++)
      {
        j = ((double)y - cy) * inc_y;
        image[n][x][y] = 0.5 * (1.0F + cos(i + (j * j)));
        
      }
    }    
  }
  image_correct_dynamic_range(image,image_height, image_width,num_of_channels);
  return (image);

}

double ***image_generator(size_t image_height, size_t image_width, int num_of_channels,double (*func)(double, double) ,double min, double max)
{
  double ***image = image_create_image(image_height, image_width,num_of_channels);
  double range = max - min;
  double inc_x = range /(double) image_height ;
  double inc_y = range /(double) image_width ;
  double cx = image_height / 2;
  double cy = image_width / 2;
  double i,j;
  for (int n = 0; n < num_of_channels; n++)
  {
    for (size_t x = 0; x < image_height; x++)
    {
      i = ((double)x - cx) * inc_x;
      i = i * i;
      for (size_t y = 0; y < image_width; y++)
      {
        j = ((double)y - cy) * inc_y;
        image[n][x][y] = (*func)(i,j);
        
      }
    }    
  }
  image_correct_dynamic_range(image,image_height, image_width,num_of_channels);
  return (image);
}
double image_circle(double x, double y)
{
  return sqrt(x * x + y * y);
}