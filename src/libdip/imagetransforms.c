
/**
 * @brief increases the difference between the highest and lowest pixel in the image.
 * @image the image.
 * @param image_height the number of rows in the image.
 * @param image_width the number of columns in the image.
 * @param num_of_channels the number of channels in the image, e.g 3 for RGB.
*/
void image_contrast_stretch(double ***image, int image_height, int image_width,int num_of_channels)
{
   double max, min,c,scaler, *table;
   table = (double *)malloc(sizeof(double) * 256);
  for (int n = 0; n < num_of_channels; n++)
  {   
    image_max_min(image[n],image_height,image_width,&max,&min);
    c = max - min;
    scaler = 255.0F/ c;
    for (int i = 0; i < 256; i++)
    {
      table[i] = scaler * ((double)i);
    }
    for (int x = 0; x < image_height; x++)
    {
      for (int y = 0; y < image_width; y++)
      {
        image[n][x][y] = table[(int)image[n][x][y]];
      }     
    }    
  }
  free(table);
  return;
}

/**
 * @brief clips the intensity values in the image, to the range [0 - 255]
 * @param image_height the number of rows in the image.
 * @param image_width the number of columns in the image.
 * @param num_of_channels the number of channels in the image, e.g 3 for RGB.
*/
void image_clip_values(double ***image, int image_height, int image_width,int num_of_channels)
{
  for (int n = 0; n < num_of_channels; n++)
  {
    for (int x = 0; x < image_height; x++)
    {
      for (int y = 0; y < image_width; y++)
      {
        if (image[n][x][y] < 0.0F)
        {
          image[n][x][y] = 0.0F;
        }else if (image[n][x][y] > 255)
        {
           image[n][x][y] = 255.0F;
        }
        
      }      
    }
  }
  return;
}

/**
 * @brief corrects the dynamic range of the image to span between [0 - 255]
 * @param image_height the number of rows in the image.
 * @param image_width the number of columns in the image.
 * @param num_of_channels the number of channels in the image, e.g 3 for RGB. 
*/
void image_correct_dynamic_range(double ***image, int image_height, int image_width,int num_of_channels)
{
  double max, min, gmax,scaler;
  for (int n = 0; n < num_of_channels; n++)
  {
    image_max_min(image[n],image_height,image_width,&max,&min);
    gmax = max - min;
    if (gmax == 0.0F)
    {
      return;
    }    
    scaler = 255.0F / gmax;
    for (int x = 0; x < image_height; x++)
    {
      for (int y = 0; y < image_width; y++)
      {
        image[n][x][y] = (image[n][x][y] - min) * scaler;
      }      
    }    
  }
  return;
}

/**
 * @brief performs histogram equalization on the image.
 * @param image_height the number of rows in the image.
 * @param image_width the number of columns in the image.
 * @param num_of_channels the number of channels in the image, e.g 3 for RGB. 
*/
void image_histogram_equalization(double ***image, int image_height, int image_width,int num_of_channels)
{
  double **stat, **s_arr;
  stat = image_get_unnormalized_histogram(image,image_height,image_width,num_of_channels);
  s_arr = malloc(sizeof(double *) * num_of_channels);
  for (int n = 0; n < num_of_channels; n++)
  {
    s_arr[n] = (double *)malloc((sizeof(double)) * 256);
    image_get_s(stat[n],s_arr[n],image_height,image_width);    
    for (int x = 0; x < image_height; x++)
    {
      for (int y = 0; y < image_width; y++)
      {
        image[n][x][y] = s_arr[n][(int)image[n][x][y]];
      }      
    }
    free(stat[n]);
    free(s_arr[n]);    
  }
  free(stat);
  free(s_arr);
  return;
}

/**
 * @brief computes the negative of the image.
 * @param image_height the number of rows in the image.
 * @param image_width the number of columns in the image.
 * @param num_of_channels the number of channels in the image, e.g 3 for RGB.
*/
void image_digital_negative(double ***image, int image_height, int image_width,int num_of_channels)
{
  double *table = (double *)malloc(sizeof(double) * 256);
  for (int i = 0; i < 256; i++)
  {
    table[i] = -(double)i + 255.0F;
  }
  for (int n = 0; n < num_of_channels; n++)
  {  
    for (int x = 0; x < image_height; x++)
    {
      for (int y = 0; y < image_width; y++)
      {
        image[n][x][y] = table[(int)image[n][x][y]];
      }      
    }        
  }
  free(table);
  return;
}

/**
 * @brief computes the log tranform of the image
 * @param image_height the number of rows in the image.
 * @param image_width the number of columns in the image.
 * @param num_of_channels the number of channels in the image, e.g 3 for RGB.
 * @param c constant to mu;tiply with the log.
*/
void image_log_tranform(double ***image, int image_height, int image_width,int num_of_channels, double c)
{
  double *table = (double *)malloc(sizeof(double) * 256);
  for (int i = 0; i < 256; i++)
  {
    table[i] = c * log10((double)i + 1.0F);
  }
  for (int n = 0; n < num_of_channels; n++)
  {  
    for (int x = 0; x < image_height; x++)
    {
      for (int y = 0; y < image_width; y++)
      {
        image[n][x][y] = table[(int)image[n][x][y]];
      }      
    }        
  }
  free(table);
  return;
}

/**
 * @brief computes the inverse log tranform of the image
 * @param image_height the number of rows in the image.
 * @param image_width the number of columns in the image.
 * @param num_of_channels the number of channels in the image, e.g 3 for RGB.
 * @param c constant to mu;tiply with the log.
*/
void image_inverselog_tranform(double ***image, int image_height, int image_width,int num_of_channels, double c)
{
  double *table = (double *)malloc(sizeof(double) * 256);
  for (int i = 0; i < 256; i++)
  {
    table[i] = c * pow(10.0F,(double)i);
  }  
  for (int n = 0; n < num_of_channels; n++)
  {  
    for (int x = 0; x < image_height; x++)
    {
      for (int y = 0; y < image_width; y++)
      {
        image[n][x][y] = table[(int)image[n][x][y]];
      }  
    }        
  }
  free(table);
  return;
}

/**
 * @brief computes the gamma tranform of the image.
 * @param image_height the number of rows in the image.
 * @param image_width the number of columns in the image.
 * @param num_of_channels the number of channels in the image, e.g 3 for RGB.
 * @param c positive constant.
 * @param gamma positive constant.
*/
void image_power_tranform(double ***image, int image_height, int image_width,int num_of_channels, double c, double gamma)
{
  double *table = (double *)malloc(sizeof(double) * 256);
  for (int i = 0; i < 256; i++)
  {
    table[i] = c * pow((double)i,gamma);
  }  
  for (int n = 0; n < num_of_channels; n++)
  {  
    for (int x = 0; x < image_height; x++)
    {
      for (int y = 0; y < image_width; y++)
      {
        image[n][x][y] = table[(int)image[n][x][y]];
      }      
    }        
  }
  free(table);
  return;
}

/**
 * @brief computes the normalized histogram equalization probabilities.
 * @param r_array array containing the unormalized frequencies for occurance of intensity values.
 * @param image_height the number of rows in the image.
 * @param image_width the number of columns in the image. 
*/
void image_get_s(double *r_array, double *s_array, int image_height,int image_width)
{
  double s, divisor;
  divisor = (double) (image_height * image_width);
  for (int k = 0; k < 256; k++)
  { 
    s = 0.0F;
    for (int j = 0; j <= k; j++)
    {
      s += r_array[j];
    }
    s = round((s * 255.0F) / divisor);
    s_array[k] = s;
  }
  return;
}

/**
 * @brief finds the highest and lowest intensity value in an image channel.
 * @param image_channel a 2D array containing the intensity values.
 * @param height number of rows in the image channel.
 * @param width number of columns in the image channel.
 * @param max memory location where the maximum intensity value is to be stored.
 * @param min memory location where the minimum intemsity value is to be stored.
*/
void image_max_min(double **image_channel, int height, int width, double *max, double *min)
{
  double mininum = image_channel[0][0];
  double maximum = image_channel[0][0];
  double pixel;
  for (int x = 0; x < height; x++)
  {
    for (int y = 0; y < width; y++)
    {
      pixel = image_channel[x][y];
      if (maximum < pixel)
      {
        maximum = pixel;
      }
      if (mininum > pixel)
      {
        mininum = pixel;
      }
      
    }
    
  }
  (*max) = maximum;
  (*min) = mininum;
  return;
}

double ***image_affine_transform(double ***image, int image_height, int image_width,int num_of_channels)
{
  typedef struct {
    int x;
    int y;
    int is_valid;
  }point;
  double ft = sin((45 / 180.0F) * 3.141592654F);
  double cf[6] = {ft,ft,-14.1421356f,-ft,ft, 0.0F};
  point **pts = malloc(sizeof(point *) * image_height);
  for (int i = 0; i < image_height; i++)
  {
    pts[i] = malloc(sizeof(point) * image_width);
  }
  int cx = image_height / 2;
  int cy = image_width / 2;
  double in_x_1;
  double in_x_2;
  double tx;
  double ty;
  for (int x = 0; x < image_height; x++)
  {
    in_x_1 = (double)(x-cx) * cf[0];
    in_x_2 = (double)(x -cx) * cf[3];

    for (int y = 0; y < image_width; y++)
    {
      pts[x][y].is_valid = 1;
      tx = in_x_1 + (double)(y -cy) * cf[1] + cf[3] + cx;
      ty = in_x_2 + (double)(y -cy) * cf[4] + cf[5] + cy;
      pts[x][y].x = (int)round(tx);
      pts[x][y].y = (int)round(ty); 
      if ((pts[x][y].x < 0) || (pts[x][y].x >= image_height) || (pts[x][y].y < 0 )|| (pts[x][y].y >= image_width))
      {
        pts[x][y].is_valid = 0;
      }
    }
    
  }
  double ***new_image = image_create_image(image_height,image_width,num_of_channels);
  for (int n = 0; n < num_of_channels; n++)
  {  
    for (int x = 0; x < image_height; x++)
    {
      for (int y = 0; y < image_width; y++)
      {
        if (pts[x][y].is_valid)
        {
          new_image[n][x][y] = image[n][pts[x][y].x][pts[x][y].y];
        }
      }      
    }        
  }
  
  for (int i = 0; i < image_height; i++)
  {
    free(pts[i]);
  }
  free(pts);
  return (new_image);
}
/**
 * @brief converts fromRGB to hsi
*/
double *** image_rgb_to_hsi (double ***rgb,int height,int width)
{
  printf("here: \n");
  double ***hsi = image_create_image(height,width,3);
  double R,G,B; 
  double pi = 3.412592654 * 180;
  double inv_pi = 1 / pi;
  double denom, numer;
  for (size_t x = 0; x < height; x++)
  {
    for (size_t y = 0; y < width; y++)
    {
      R = rgb[0][x][y];
      G = rgb[1][x][y];
      B = rgb[2][x][y];
      hsi[2][x][y] = (R + G + B) / 3.0;
      if (((int)R == 0) && ((int)G == 0) && ((int)B == 0))
      {
        hsi[0][x][y] = 0;
        hsi[1][x][y] = 0;

      }else
      {
          denom = sqrt((R - G) * (R - G) + (R - B) * (G - B));
          numer = 0.5 * ((R - G) + (R - B));
          if ((int)B <= (int)G)
          {
            hsi[0][x][y] = acos(numer / denom) * inv_pi;
            if (R < B)
            {
              hsi[1][x][y] = 1.0 - 3.0 * R / (R + G + B);
            }else
            {
              hsi[1][x][y] = 1.0 - 3.0 * B / (R + G + B);
            }
          } else
          {
            hsi[0][x][y] = 360.0 - acos(numer / denom) * inv_pi; 
            if (R < G)
            {
              hsi[1][x][y] = 1.0 - 3.0 * R / (R + G + B);
            }else
            {
              hsi[1][x][y] = 1.0 - 3.0 * G / (R + G + B);
            }
         }                      
      }     
    }
  } 
  return (hsi);
}
/**
 * @brief converts from hsi to RGB
*/
double *** image_hsi_to_rgb (double ***hsi,int height,int width)
{
  double ***rgb = image_create_image(height,width,3);
  double R , G , B, H , S, I;
  double pi = 3.412592654 / 180;
  for (size_t x = 0; x < height; x++)
  {
    for (size_t y = 0; y < width; y++)
    {
      H = hsi[0][x][y];
      S = hsi[1][x][y];
      I = hsi[2][x][y];
      if (H < 120.0)
      {
        B = I * ( 1 - S );
        R = I * (1 + S * cos ( H * pi ) / cos((60 - H) * pi));
        G = 3.0 * I - ( R + B);
      }else if( H < 240.0)
      {
        H -= 120.0;
        R = I * ( 1 - S );
        G = I * (1 + S * cos ( H * pi ) / cos((60 - H) * pi));
        B = 3.0 * I - ( R + G);
      }else
      {
        H -= 240.0;
        G = I * ( 1 - S );
        B = I * (1 + S * cos ( H * pi ) / cos((60 - H) * pi));
        R = 3.0 * I - ( G + B);
      } 
      rgb[0][x][y] = R;
      rgb[1][x][y] = G;
      rgb[2][x][y] = B;
    }
  }
  return (rgb);
}