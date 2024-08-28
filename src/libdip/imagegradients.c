double*** image_laplacian(double ***image,size_t image_height, size_t image_width, int num_of_channels)
{
  double kernel[3][3] = {{1,1,1},{1,-8,1},{1,1,1}};
  //double kernel[3][3] = {{0,1,0},{1,-4,1},{0,1,0}};
  double ***new_image = image_create_image(image_height,image_width,num_of_channels);
  int a,b;
  a = b = 1;
  int i, j;
  double gxy;
  for (size_t n = 0; n < num_of_channels; n++)
  {
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        gxy = 0.0F;
        for (int s = -a; s <= a; s++)
        {
          i = x - s;
          //perform mirror padding.
          if (i < 0)
          {
            i = -i; // reflect about the x axis.
          } else if (i >= image_height)
          {
            i = x + s;
          }
          for (int t = -b; t <= b; t++)
          {
            j = y - t;
            if (j < 0)
            {
              j = -j; //reflect about the x axis.
            } else if (j >= image_width)
            {
              j = y + t;
            }
            gxy += kernel[s + a][ t + b] * image[n][i][j];        
          } 
        }
        new_image[n][x][y] = gxy;
      }      
    }
  }
  return (new_image);
}

void image_sharpen_laplacian(double ***image,size_t image_height, size_t image_width, int num_of_channels, double c)
{
  double ***new_image = image_laplacian(image,image_height,image_width,num_of_channels);
  for (size_t n = 0; n < num_of_channels; n++)
  {
    //convolve along rows.
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        image[n][x][y] = image[n][x][y] + c * new_image[n][x][y];
      }      
    }
  }
  image_free_image(new_image,image_height,num_of_channels);
  image_clip_values(image,image_height,image_width,num_of_channels);
  return;
}

double ***image_sobel_gradient(double ***image,size_t image_height, size_t image_width, int num_of_channels)
{
  double gx_kernel[3][3] = {{-1,-2,-1},{0,0,0},{1,2,1}};
  double gy_kernel[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
  double ***new_image = image_create_image(image_height,image_width,num_of_channels);
  int a,b;
  a = b = 1;
  int i, j;
  double gx,gy;
  for (size_t n = 0; n < num_of_channels; n++)
  {
    //convolve along rows.
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        gx = 0.0F;
        gy = 0.0F;
        // finding gx.
        for (int s = -a; s <= a; s++)
        {
          i = x - s;
          //perform mirror padding.
          if (i < 0)
          {
            i = -i; // reflect about the x axis.
          } else if (i >= image_height)
          {
            i = x + s;
          }
          for (int t = -b; t <= b; t++)
          {
            j = y - t;
            if (j < 0)
            {
              j = -j; //reflect about the x axis.
            } else if (j >= image_width)
            {
              j = y + t;
            }
            gx += gx_kernel[s + a][ t + b] * image[n][i][j];        
          } 
        }
        // finding gy.
        for (int s = -a; s <= a; s++)
        {
          i = x - s;
          //perform mirror padding.
          if (i < 0)
          {
            i = -i; // reflect about the x axis.
          } else if (i >= image_height)
          {
            i = x + s;
          }
          for (int t = -b; t <= b; t++)
          {
            j = y - t;
            if (j < 0)
            {
              j = -j; //reflect about the x axis.
            } else if (j >= image_width)
            {
              j = y + t;
            }
            gy += gy_kernel[s + a][ t + b] * image[n][i][j];        
          } 
        }
        new_image[n][x][y] = sqrt(gx * gx + gy * gy);
      }      
    }
  }
  return (new_image);
}