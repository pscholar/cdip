/**
 * @brief performs spatial correlation.
 * @param image_channel  array containing the pixels of an image channel
 * @param kernel  array containing the kernel coefficients, size of kernel m x n.
 * @param image_height number of rows in the image channel.
 * @param image_width number of columns in the image channel.
 * @param cx displacement along the x-axis.
 * @param cy displacement along the y-axis.
 * @param a a = (m - 1) / 2.
 * @param b b = (n - 1) / 2.
 * 
*/
double image_spatial_gxy_corr(double **image_channel,double **kernel,size_t image_height, size_t image_width,int cx, int cy, int a, int b )
{
  double gxy = 0.0F;
  int x;
  int y;
  int is_x_valid;
  int is_y_valid;

  for (int s = -a; s <= a; s++)
  {
    x = cx + s;
    is_x_valid = (x >= 0) && (x < image_height);
    if(is_x_valid)
    {
      for (int t = -b; t <= b; t++)
      {
        y = cy + t;
        is_y_valid = (y >= 0) && (y < image_width);
        if (is_y_valid)
        {
          gxy += kernel[s + a][ t + b] * image_channel[x][y];
        }      
         
      }
    }
    
    
  }
  return (gxy);
}
/**
 * @brief performs spatial convolution.
 * @param image_channel  array containing the pixels of an image channel
 * @param kernel  array containing the kernel coefficients, size of kernel m x n.
 * @param image_height number of rows in the image channel.
 * @param image_width: number of columns in the image channel.
 * @param cx displacement along the x-axis.
 * @param cy displacement along the y-axis.
 * @param a a = (m - 1) / 2.
 * @param b b = (n - 1) / 2.
 * 
*/
double image_spatial_gxy_conv(double **image_channel,double **kernel,size_t image_height, size_t image_width,int cx, int cy, int a, int b )
{
  double gxy = 0.0F;
  int x;
  int y;
  for (int s = -a; s <= a; s++)
  {
    x = cx - s;
    //perform mirror padding.
    if (x < 0)
    {
      x = -x; // reflect about the x axis.
    } else if (x >= image_height)
    {
      x = cx + s;
    }
    for (int t = -b; t <= b; t++)
    {
      y = cy - t;
      if (y < 0)
      {
        y = -y; //reflect about the x axis.
      } else if (y >= image_width)
      {
        y = cy + t;
      }
      gxy += kernel[s + a][ t + b] * image_channel[x][y];        
    } 
  }
  return (gxy);
}
/**
 * @brief performs linear spatial filtering
 * @param image array containing image data.
 * @param kernel array containing kernel coefficients.
 * @param image_height number of rows in the image.
 * @param image_width number of columns in the image.
 * @param num_of_channels number of channels in the image.
 * @param m number of rows in the kernel.
 * @param n number of columns in the kernel.
*/
double ***image_spatial_filter(double ***image,double **kernel,size_t image_height, size_t image_width, int num_of_channels, int m , int n)
{
  double ***new_image = image_create_image(image_height,image_width, num_of_channels);
  int a = (m - 1) / 2;
  int b =(n - 1) / 2;
  for (size_t n = 0; n < num_of_channels; n++)
  {
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        new_image[n][x][y]  = image_spatial_gxy_corr(image[n],kernel,image_height,image_width,x,y,a,b);       
      }      
    }
  }
  //image_free_image(image,image_height,num_of_channels);
  return (new_image);
}
double ***image_spatial_smooth(double ***image,size_t image_height, size_t image_width, int num_of_channels, int m , int n,double k)
{
  double **kernel = image_create_box_kernel_un(m,n);
  //double **kernel = image_create_gaussian_kernel_un(m,n,k);
  double ***new_image = image_spatial_filter(image,kernel,image_height,image_width,num_of_channels,m,n);
  image_free_kernel(kernel,m);
  return (new_image);
}


/*Implementation for a separable kernel*/
double ***image_spatial_filter_sep(double ***image,double *kernel_row,double *kernel_col,size_t image_height, size_t image_width, int num_of_channels, int m , int n)
{
  double ***new_image = image_create_image(image_height,image_width, num_of_channels);
  int a = (m - 1) / 2;
  int b =(n - 1) / 2;
  double gxy;
  int i;
  int j;
  for (size_t n = 0; n < num_of_channels; n++)
  {
    //convolve along rows.
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        gxy = 0.0F; 
        for (int s = -a; s <= a; s++)
        {
          i = y - s;
          // perform mirror padding.
          if (i >= image_width)
          {
            i = y + s;
          }else if (i < 0 )
          {
            i = -i;
          }  
          gxy += kernel_row[s + a] * image[n][x][i];
        }
        new_image[n][x][y]  = gxy; 
      }      
    }
  }
  //image_free_image(image,image_height,num_of_channels);
  double ***new_image_2 = image_create_image(image_height,image_width, num_of_channels);
  //convolve along columns.
  for (size_t n = 0; n < num_of_channels; n++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      for (size_t x = 0; x < image_height; x++)
      {
        gxy = 0.0F; 
        for (int t = -b; t <= b; t++)
        {
          j = x - t;
          // perform mirror padding.
          if (j >= image_height)
          {
            j = x + t;
          }else if (j < 0 )
          {
            j = -j;
          }  
          gxy += kernel_col[t + b] * new_image[n][j][y];       
        }
        new_image_2[n][x][y]  = gxy; 
      }      
    }
  }
  image_free_image(new_image,image_height,num_of_channels);
  return (new_image_2);
}

double ***image_spatial_smooth_sep(double ***image,size_t image_height, size_t image_width, int num_of_channels, int m , int n,double k)
{
  //double *kernel = image_create_box_kernel_sep(m);
  double *kernel = image_create_gaussian_kernel_sep(m,1.0);
  double ***new_image = image_spatial_filter_sep(image,kernel,kernel,image_height,image_width,num_of_channels,m,n);
  return (new_image);
}

double ***image_unsharpmask(double ***image,size_t image_height, size_t image_width, int num_of_channels, int m , int n,double k, double C)
{
  double ***new_image = image_spatial_smooth_sep(image,image_height,image_width,num_of_channels,m,n,k);
  for (size_t n = 0; n < num_of_channels; n++)
  {
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        new_image[n][x][y] = image[n][x][y] + C * (image[n][x][y] - new_image[n][x][y]);
      }      
    }
  }
  image_clip_values(new_image, image_height,image_width,num_of_channels);
  return (new_image);
}

double ***image_order_static_filter(double ***image,size_t image_height, size_t image_width, int num_of_channels, int m,int t)
{
  double *arr = (double *)malloc(sizeof(double) * m * t);
  double ***new_image = image_create_image(image_height,image_width,num_of_channels);
  for (int n = 0; n < num_of_channels; n++)
  {
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        new_image[n][x][y] =  image_order_values(image[n],image_height,image_width,x,y,arr,m,t);
      }    
    }
  }
  free(arr);
  image_free_image(image,image_height,num_of_channels);
  return (new_image);
}

double image_order_values(double **image_channel,size_t height, size_t width,int cx, int cy,double *arr, int m, int n)
{
  int a = (m - 1) / 2;
  int b =(n - 1) / 2;
  int x;
  int y;
  int i = 0;
  for (int s = -a; s <= a; s++)
  {
    x = cx - s;
    if (x < 0)
    {
      x = -x;
    } else if (x >= height)
    {
      x = cx + s;
    }
    for (int t = -b; t <= b; t++)
    {
      y = cy - t;
      if (y < 0)
      {
        y = -y;
      } else if (y >= width)
      {
        y = cy + t;
      }
      arr[i] = image_channel[x][y];
      i++;
    }  
  }
  support_quicksort(arr,0, i - 1);
  int index = (i / 2);
  return (arr[index]);
}