/**
 * @brief performs spatial convolution with an unseparable kernel
*/
double** image_seg_spatial_convolution(double **image, size_t image_height, size_t image_width,int m, int n, double (*kernel)[n])
{
  //using replicative padding
  double** nimage = image_create_image_channel(image_height,image_width);
  int tx,ty, a,b;
  a = m / 2;
  b = n / 2;
  double weighted_value;
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      weighted_value = 0.0F;
      for (int s = -a; s <= a; s++)
      {
        tx = x + s;
        
        if (tx >= image_height)
        {
          tx = x - s;
        }else if (tx < 0)
        {
          tx = -tx;
        }
        for (int t = -b; t <= b; t++)
        {
          ty = y + t;
          if (ty >= image_width)
          {
            ty = y - t;
          }else if(ty < 0){
            ty = -ty;
          }
          weighted_value += image[tx][ty] * kernel[s + a][t + b];
        }
      } 
      nimage[x][y] = weighted_value;    
    }    
  }
  return (nimage);  
}

/**
 * @brief detects a point using laplacian
*/
double** image_seg_laplacian_point(double **image, size_t image_height, size_t image_width,int m, int n, double (*kernel)[n], double fraction)
{
  double** lap_image = image_seg_spatial_convolution(image,image_height,image_width,m,n,kernel);
  double min, max;
  image_max_min(lap_image,image_height,image_width,&max,&min);
  printf("max: %d : min %d\n",(int)max,(int)min);
  min = min < 0 ? -min : min;
  max = max < 0 ? -max : max;
  max = max > min ? max : min;
  double threshold = fraction * max;
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      lap_image[x][y] = lap_image[x][y] < 0 ? -lap_image[x][y] : lap_image[x][y];
      lap_image[x][y] = lap_image[x][y] > threshold ? 255.0F : 0.0F;
    }    
  }
  return (lap_image);
}
/**
 * @brief detects lines in an image
*/
double** image_seg_laplacian_line_detection(double **image, size_t image_height, size_t image_width,int m, int n, double (*kernel)[n], double fraction)
{
  double** lap_image = image_seg_spatial_convolution(image,image_height,image_width,m,n,kernel);
  double min, max;
  image_max_min(lap_image,image_height,image_width,&max,&min);
  double threshold = max * fraction;
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      lap_image[x][y] = lap_image[x][y] > threshold? 255.0F : 0.0F;
    }    
  }
  return (lap_image);
}

double** image_seg_g(double **image, size_t image_height, size_t image_width,int m, int n, double (*kernel)[n])
{
  double **g = image_seg_spatial_convolution(image,image_height,image_width,m,n,kernel);
  return (g);
}
/**
 * calculates the gradient image
*/
double** image_seg_gradient_image(double **image, size_t image_height, size_t image_width,int m, int n, double (*kernel_x)[n],double (*kernel_y)[n])
{
  double **gx = image_seg_spatial_convolution(image,image_height,image_width,m,n,kernel_x);
  double **gy = image_seg_spatial_convolution(image,image_height,image_width,m,n,kernel_y);
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      gx[x][y] = gx[x][y] < 0.0F? -gx[x][y] : gx[x][y];
      gy[x][y] = gy[x][y] < 0.0F? -gy[x][y] : gy[x][y];
      gx[x][y] = gx[x][y] + gy[x][y];
    }    
  }
  image_free_kernel(gy,image_height);
  return (gx);
}
/**
 * @brief calculates the gradient direction
*/
double** image_seg_gradient_direction(double **image, size_t image_height, size_t image_width,int m, int n, double (*kernel_x)[n],double (*kernel_y)[n])
{
  double **gx = image_seg_spatial_convolution(image,image_height,image_width,m,n,kernel_x);
  double **gy = image_seg_spatial_convolution(image,image_height,image_width,m,n,kernel_y);
  double scaler = 180.0 / 3.141592654;
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      gx[x][y] = scaler * atan2(gy[x][y] , gx[x][y] );
    }    
  }
  image_free_kernel(gy,image_height);
  return (gx);
}

/**
 * @brief computes a samples of a gaussian through it's center
 * 
*/
double *image_seg_create_gaussian_kernel_(double std, int *m)
{
  int max = (int)ceil(6 * std);
  
  max = abs(max);
  if(!(max % 2))
  {
    max = max + 1;
  }
  *(m) = max;
  int c = max / 2;
  double* gaus = (double *)malloc(sizeof(double) * max);
  double mul = -(2.0) * std * std;
  double sum = 0.0;
  for (int t = -c; t <= c; t++)
  {
    gaus[t + c] = exp( ((double)( t * t) )/ (mul));
    sum += gaus[t + c];
  }
  printf("m = %d\n",*m);
  for (int i = 0; i < max; i++)
  {
    gaus[i] /= sum;
  }
  return (gaus);
}
/**
 * @brief performs spatial convolution of a separable kernel
*/

double **image__seg_spatial_convolution_sep(double **image,double *kernel_row,double *kernel_col,size_t image_height, size_t image_width, int m , int n)
{
  int a,b,i,j;
  double gxy, **new_image_rows = image_create_image_channel(image_height,image_width);
  a = (m - 1) / 2;
  b = (n - 1) / 2;
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
        gxy += kernel_row[s + a] * image[x][i];
      }
      new_image_rows[x][y]  = gxy; 
    }      
  }
  double **new_image_cols = image_create_image_channel(image_height,image_width);
  //convolve along columns.
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
        gxy += kernel_col[t + b] * new_image_rows[j][y];       
      }
      new_image_cols[x][y]  = gxy; 
    }      
  }
  
  image_free_kernel(new_image_rows,image_height);
  return (new_image_cols);
}

/**
 * @brief marks the non zero crossing in a laplacian image
*/
double** image_seg_laplacian_zero_crossing(double **image, size_t image_height, size_t image_width,double fraction)
{
  int topx,bottomx,lefty,righty;
  double min, max;
  image_max_min(image,image_height,image_width,&max,&min);
  double threshold = max * fraction;
  printf("max: %f : min: %f :threhold: %f\n",max,min,threshold);
  double diff;
  double **zero_cross_image = image_create_image_channel(image_height,image_width);
  //convolve along rows.
  for (size_t x = 0; x < image_height; x++)
  {
    topx = x - 1;
    bottomx = x + 1;
    if ((topx < 0)|| (bottomx >= image_height))
    {
      continue;
    }
    
    for (size_t y = 0; y < image_width; y++)
    { 
      lefty = y - 1;
      righty = y + 1;
      if ((lefty < 0 ) || (righty >= image_width))
      {
        continue;
      }
      if (((image[topx][y] < 0.0 && image[bottomx][y] > 0.0)) || ((image[topx][y] > 0.0 && image[bottomx][y] < 0.0)))
      {
        diff = image[topx][y] + image[bottomx][y];
        diff = diff > 0.0 ? diff: -diff;
        if (diff > threshold)
        {
          zero_cross_image[x][y] = 255.0;
          continue;
        }        
      }
      if (((image[x][lefty] < 0.0 && image[x][righty] > 0.0)) || ((image[x][lefty] > 0.0 && image[x][righty] < 0.0)))
      {
        diff = image[x][lefty] + image[x][righty];
        diff = diff > 0.0 ? diff: -diff;
        if (diff > threshold)
        {
          zero_cross_image[x][y] = 255.0;
          continue;
        }        
      }
      if (((image[topx][lefty] < 0.0 && image[bottomx][righty] > 0.0)) || ((image[topx][lefty] > 0.0 && image[bottomx][righty] < 0.0)))
      {
        diff = image[topx][lefty] + image[bottomx][righty];
        diff = diff > 0.0 ? diff: -diff;
        if (diff > threshold)
        {
          zero_cross_image[x][y] = 255.0;
          continue;
        }        
      }
      if (((image[bottomx][lefty] < 0.0 && image[topx][righty] > 0.0)) || (image[bottomx][lefty] > 0.0 && image[topx][righty] < 0.0))
      {
        diff = image[bottomx][lefty] + image[topx][righty];
        diff = diff > 0.0 ? diff: -diff;
        if (diff > threshold)
        {
          zero_cross_image[x][y] = 255.0;
          continue;
        }        
      }
      
    }      
  }
  return (zero_cross_image);
}

/**
 * @brief implements Marr_Hildereth edge detector
*/
double** image_seg_marr_hildereth_edge(double **image, size_t image_height, size_t image_width,int m, int n, double (*lap_kernel)[n],double std, double fraction)
{
  //1. smooth the image with a gaussian kernel
  int size;
  double *gauss_kernel = image_seg_create_gaussian_kernel_(std, &size);
  double **smoothed_image = image__seg_spatial_convolution_sep(image,gauss_kernel,gauss_kernel,image_height,image_width,size,size);
  free(gauss_kernel);
  //2. compute the laplacian of the smoothed image
  double **laplacian_image = image_seg_spatial_convolution(smoothed_image,image_height,image_width,m,n,lap_kernel);  
  image_free_kernel(smoothed_image,image_height);
  //3. find the zero crossing of the image.
  double **edge_image = image_seg_laplacian_zero_crossing(laplacian_image,image_height,image_width,fraction);
  image_free_kernel(laplacian_image,image_height);
  return (edge_image);
}

/**
 * @brief implementation of the canny edge detector
*/

double** image_seg_canny_edge(double **image, size_t image_height, size_t image_width,int m, int n, double (*kernel_x)[n], double (*kernel_y)[n],double std ,double t_low, double t_high)
{
  //1. smooth the image
  int size;
  double *gauss_kernel = image_seg_create_gaussian_kernel_(std, &size);
  double **smoothed_image = image__seg_spatial_convolution_sep(image,gauss_kernel,gauss_kernel,image_height,image_width,size,size);
  free(gauss_kernel);
  //2.compute the gradient image and gradient direction
  double **gx = image_seg_spatial_convolution(smoothed_image,image_height,image_width,m,n,kernel_x);
  double **gy = image_seg_spatial_convolution(smoothed_image,image_height,image_width,m,n,kernel_y);
  image_free_kernel(smoothed_image,image_height);
  double **mag = gx;
  double **angle = image_create_image_channel(image_height, image_width);
  double degree= 180.0 / 3.141592654; 
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      angle[x][y] = degree * atan2(gy[x][y] , gx[x][y]);
      gx[x][y] = gx[x][y] < 0.0? -gx[x][y] : gx[x][y];
      gy[x][y] = gy[x][y] < 0.0? -gy[x][y] : gy[x][y];
      mag[x][y] = gx[x][y] + gy[x][y];
    }    
  }
  image_free_kernel(gy,image_height);
  //3. perfom nonlocal maxima suppression
  double **gn = image_create_image_channel(image_height, image_width); 
  int direction;
  int txtop,tx_bottom, tyleft,tyright;
  for (size_t x = 0; x < image_height; x++)
  {
    txtop = x - 1;
    tx_bottom = x + 1;
    if (txtop < 0)
    {
      txtop = x + 1;
    }
    if (tx_bottom >= image_height)
    {
      tx_bottom = x - 1;
    }   
    for (size_t y = 0; y < image_width; y++)
    {
      tyleft = y - 1;
      tyright = y + 1;
      if (tyleft < 0)
      {
        tyleft = y + 1;
      }
      if (tyright >= image_width)
      {
        tyright = y - 1;
      }
      direction = (int)image_seg_map_angle(angle[x][y]);
      //printf("txtop:%d txbottom : %d tyleft: %d tyright: %d\n",txtop,tx_bottom,tyleft,tyright);
      switch (direction)
      {
        case 0:
          if ((mag[x][y] >= mag[txtop][y]) && (mag[x][y] >= mag[tx_bottom][y])) 
          {
            gn[x][y] = mag[x][y];
          }       
          break;
        case -45:
          if ((mag[x][y] >= mag[txtop][tyleft]) && (mag[x][y] >= mag[tx_bottom][tyright])) 
          {
            gn[x][y] = mag[x][y];
          } 
          break;
        case 90:
          if ((mag[x][y] >= mag[x][tyleft]) && (mag[x][y] >= mag[x][tyright])) 
            {
              gn[x][y] = mag[x][y];
            } 
          break;
        case 45:
          if ((mag[x][y] >= mag[txtop][tyright]) && (mag[x][y] >= mag[tx_bottom][tyleft])) 
          {
            gn[x][y] = mag[x][y];
          } 
          break;
      }
    }    
  }
  image_free_kernel(mag,image_height);
  image_free_kernel(angle,image_height); 
  //4. link edge points using hysteresis thresholding
  double **canny = image_seg_hysteresis_threshold(gn,image_height,image_width,t_low,t_high);
  return (canny);
}
/**
 * @brief maps an angle to a particular direction
*/
double image_seg_map_angle(double angle)
{
  double m_angle ;
  m_angle = angle < 0.0? angle + 180.0 : angle;
  if ((m_angle >= 0 && m_angle < 22.5) || (m_angle >= 157.5))
  {
    m_angle = 0.0; // horizontal
  }else if ((m_angle >= 22.5) && (m_angle < 67.5))
  {
    m_angle = -45.0;
  }else if ((m_angle >= 67.5) && (m_angle < 112.5))
  {
    m_angle = 90.0; // vertical
  }else if ((m_angle >= 112.5) && (m_angle < 157.5))
  {
    m_angle = +45.0;
  } else
  {
    m_angle = 0.0;
  }
  return (m_angle);
}
/**
 * @brief performs hysteresis thresholding of an image.
*/
double** image_seg_hysteresis_threshold(double **image, size_t image_height, size_t image_width, double t_low, double t_high)
{
  double **glow = image_create_image_channel(image_height,image_width);
  double **ghigh = image_create_image_channel(image_height,image_width);
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      if (image[x][y] >= t_high)
      {
        ghigh[x][y]= image[x][y];        
      } else if (image[x][y] >= t_low)
        {
         glow[x][y]= image[x][y];
        }     
    }    
  }
  image_free_kernel(image,image_height);
  double **marker = image_create_image_channel(image_height,image_width);
  int tx,ty;
  int a =3/2;
  int b = 3/2;
   for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      if (ghigh[x][y] <= 0.0)
      {
        continue;
      }     
      for (int s = -a; s <= a; s++)
      {
        tx = x + s;        
        if (tx >= image_height)
        {
          tx = x - s;
        }else if (tx < 0)
        {
          tx = -tx;
        }
        for (int t = -b; t <= b; t++)
        {
          ty = y + t;
          if (ty >= image_width)
          {
            ty = y - t;
          }else if(ty < 0){
            ty = -ty;
          }
          if (glow[tx][ty] > 0.0)
          {
            marker[tx][ty] = 600.0;
          }         
        }
      }     
    }    
  }
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      if ((((int)marker[x][y]) == 600) && (ghigh[x][y] <= 0 ))
      {
        ghigh[x][y]= glow[x][y];        
      }     
    }    

  }
  image_free_kernel(glow,image_height);
  image_free_kernel(marker,image_height);
  return (ghigh);
}

/**
 * @brief performs edge linking using local processing
*/
double** image_seg_local_edge_linking(double **image, size_t image_height, size_t image_width,int m, int n, double (*kernel_x)[n], double (*kernel_y)[n],double tm, double ta, double r,int gap_length)
{
  //1. compute the gradient image and gradient direction
  double **gx = image_seg_spatial_convolution(image,image_height,image_width,m,n,kernel_x);
  double **gy = image_seg_spatial_convolution(image,image_height,image_width,m,n,kernel_y);
  double **mag = gx;
  double **angle = image_create_image_channel(image_height, image_width);
  double degree= 180.0 / 3.141592654; 
  int upper_angle_limit = (int)(ta + r);
  int lower_angle_limit = (int)(ta - r);
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      angle[x][y] = degree * atan2(gy[x][y] , gx[x][y]);
      angle[x][y] = angle[x][y] < 0.0? angle[x][y] + 180.0 : angle[x][y];   
      gx[x][y] = gx[x][y] < 0.0? -gx[x][y] : gx[x][y];
      gy[x][y] = gy[x][y] < 0.0? -gy[x][y] : gy[x][y];
      mag[x][y] = gx[x][y] + gy[x][y];
    }    
  }
  //image_free_kernel(angle,image_height);
  image_free_kernel(gy,image_height);
  double max , min;
  image_max_min(mag,image_height,image_width,&max,&min);
  double** bimage = image_create_image_channel(image_height,image_width);  
  double Tm = max * tm;
  printf("TM: %f : max: %f : min : %f\n",Tm,max,min);
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width; y++)
    {
      if ((mag[x][y] > Tm) && ((((int)angle[x][y]) <= upper_angle_limit) || (((int)angle[x][y]) >= lower_angle_limit)))
      {
        bimage[x][y] = 1;
      }      
    }    
  }
  image_free_kernel(mag,image_height);
  image_free_kernel(angle,image_height);
  return (bimage);
  int start_of_zero;
  int end_of_zero;
  int n_of_zeros;
  double **vert = image_create_image_channel(image_height,image_width);
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width;y++)
    {     
        vert[x][y] = bimage[x][y];
    }
    
  }
  // scan along rows
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width;)
    {
        if ((int)bimage[x][y] == 1)
        {
          //move until you incur a zero
          y++;
          while ((y < image_width) && (int)bimage[x][y] == 1)
          {
            y++;
          }
          if ((y < image_width))
          {
            // a gap was found, count the gaps
            start_of_zero = y;
            y++;
            while ((y < image_width) && (int)bimage[x][y] == 0)
            {
              y++;
            }
            // either we have reached the end or a 1 has been found.
            if (y < image_width)
            {
              end_of_zero = y;
              n_of_zeros = end_of_zero - start_of_zero;
              if (n_of_zeros <= gap_length)
              {
                //fillup the zeros
                for (int i = start_of_zero; i < end_of_zero; i++)
                {
                  bimage[x][i] = 1.0;
                }                
              }
            }           
          }         
        }else
        {
          y++;
        }  
    }    
  }
  // scan along columns
  for (int y = 0; y < image_width; y++)
  {
    for (int x = 0; x < image_height;)
    {
        if ((int)vert[x][y] == 1)
        {
          //move until you incur a zero
          x++;
          while ((x < image_height) && (int)vert[x][y] == 1)
          {
            x++;
          }
          if ((x < image_height))
          {
            // a gap was found, count the gaps
            start_of_zero = x;
            x++;
            while ((x < image_height) && (int)vert[x][y] == 0)
            {
              x++;
            }
            // either we have reached the end or a 1 has been found.
            if (x < image_height)
            {
              end_of_zero = x;
              n_of_zeros = end_of_zero - start_of_zero;
              if (n_of_zeros <= gap_length)
              {
                //fillup the zeros
                for (int i = start_of_zero; i < end_of_zero; i++)
                {
                  vert[i][y] = 1.0;
                }                
              }
            }           
          }         
        }else
        {
          x++;
        }  
    }    
  } 
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width;y++)
    {
      if ((int)vert[x][y] == 1)
      {
        bimage[x][y] = 1.0;
      }
    }
  }
  image_free_kernel(vert,image_height);
  return (bimage);
}
/**
 * @brief performs global thresholding on an image
*/
double** image_seg_global_thres(double **pt_image, int height, int width, double delta, int max_inter)
{
  double prev_thres, curr_thres, diff, mean_one, mean_two,val[256] = {0},average_intes = 0,pb = 1.0 / ((double)(height * width));
  double prob_one, prob_two;
  int num_inter = 0;
  for (int x = 0; x < height; x++)
  {
    for (int y = 0; y < width;y++)
    {
      val[(int)pt_image[x][y]] += pb;
    }
  }
  for (size_t i = 0; i < 256; i++)
  {
    average_intes += ((double)i) * val[i]; 
  }

  prev_thres = average_intes;
  int k = (int)round(prev_thres);
  do
  {
    mean_one = 0; mean_two = 0; prob_one = 0; prob_two = 0;
    for (size_t i = 0; i < k; i++)
    {
      mean_one += ((double)i) * val[i];
      prob_one += val[i];
    }
    prob_two = 1 - prob_one;
    for (size_t j = k; j < 256; j++)
    {
      mean_two += ((double)j) * val[j];
    }
    if (prob_one > 0)
    {
      mean_one /= prob_one;
    }
    if (prob_two > 0)
    {
      mean_two /= prob_two;
    }    
    curr_thres = (mean_one + mean_two) / 2.0;
    diff = (curr_thres - prev_thres) > 0 ? (curr_thres - prev_thres) : -(curr_thres - prev_thres);
    prev_thres = curr_thres; 
    k = (int)round(prev_thres); 
    num_inter++;
  } while ((diff > delta) && (num_inter < max_inter));  
  printf("Iterations: %3d : Threshold: %f : Average: %f\n",num_inter,prev_thres,average_intes);
  double **pt_new_image = image_create_image_channel(height,width);
  for (int x = 0; x < height; x++)
  {
    for (int y = 0; y < width;y++)
    {
      if (pt_image[x][y] < k)
      {
        pt_new_image[x][y] = 255.0;
      }else
      {
        pt_new_image[x][y] = 127.0;
      }     
    }
  }
  return (pt_new_image);
}
/**
 * @brief implements otsu's thresholding method
*/
double** image_seg_otsu_thres(double **pt_image, int height, int width)
{
  double histogram[256] = {0};
  double global_av = 0;
  double prob_one, prob_two;
  double pb = 1.0 / ((double)(height * width));
  for (int x = 0; x < height; x++)
  {
    for (int y = 0; y < width;y++)
    {
      histogram[(int)pt_image[x][y]] += pb;
    }
  }
  for (size_t i = 0; i < 256; i++)
  {
    global_av += ((double)i) * histogram[i]; 
  }
  double in_btw_var[256] = {0};
  double cum_mean[256] = {0};
  double cum_prob[256] = {0};
  double threshold;
  cum_mean[0] = 0;
  cum_prob[0] = histogram[0];
  prob_two = 1 - cum_prob[0];
  if (cum_prob[0] > 0 && prob_two > 0)
  {
    in_btw_var[0] = ((global_av * cum_prob[0] - cum_mean[0]) * (global_av * cum_prob[0] - cum_mean[0])) / (cum_prob[0] * (prob_two));
  }
  double max = in_btw_var[0];
  threshold = 0; 
  for (size_t k = 1; k < 256; k++)
  {
    cum_mean[k] = cum_mean[k - 1] + ((double)k )* histogram[k];
    cum_prob[k] = cum_prob[k - 1] + histogram[k];
    prob_two = 1 - cum_prob[k];
    if (cum_prob[k] > 0 && prob_two > 0)
    {
      in_btw_var[k] = (global_av * cum_prob[k] - cum_mean[k]) * (global_av * cum_prob[k] - cum_mean[k]) / (cum_prob[k] * prob_two);
      if (in_btw_var[k] > max)
      {
        max = in_btw_var[k];
        threshold = k;
      }else if (in_btw_var[k] == max)
      {
        threshold = (threshold + k) / 2;
      }
            
    }
  }
  printf("\n");
  printf("Threshold: %f : Max: %f \n",threshold, max);
  double **pt_new_image = image_create_image_channel(height,width);
  for (int x = 0; x < height; x++)
  {
    for (int y = 0; y < width;y++)
    {
      if (pt_image[x][y] < threshold)
      {
        pt_new_image[x][y] = 255.0;
      }else
      {
        pt_new_image[x][y] = 0;
      }     
    }
  }
  return (pt_new_image);
}

/**
 * @brief implements otsu's thresholding method
*/
double** image_seg_otsu_dualthres(double **pt_image, int height, int width)
{
  double histogram[256] = {0};
  double global_av = 0;
  double pb = 1.0 / ((double)(height * width));
  for (int x = 0; x < height; x++)
  {
    for (int y = 0; y < width;y++)
    {
      histogram[(int)pt_image[x][y]] += pb;
    }
  }
  for (size_t i = 0; i < 256; i++)
  {
    global_av += ((double)i) * histogram[i]; 
  }
  double in_btw_var[256][256] = {{0}};
  double cum_mean[256] = {0};
  double cum_prob[256] = {0};
  double threshold_one = 0, threshold_two = 0;
  double max = 0;
  cum_mean[0] = 0;
  cum_prob[0] = histogram[0]; 
  for (size_t i = 1; i < 256; i++)
  {
    cum_prob[i] = cum_prob[i - 1] + histogram[i];
    cum_mean[i] = cum_mean[i - 1] + ((double)i) * histogram[i];
  }
  double mean_one, mean_two, mean_three, prob_one, prob_two , prob_three;
  for (size_t k = 1; k < 253; k++)
  {
    prob_one = cum_prob[k];
    if (prob_one > 0)
    {
      mean_one = cum_mean[k] / prob_one;
    }else
    {
      mean_one = 0;
    }
    for (size_t i = k + 1; i < 254; i++)
    {
      prob_two = cum_prob[i] - prob_one;
      if (prob_two > 0)
      {
        mean_two = (cum_mean[i] - cum_mean[k]) / prob_two;
      }else
      {
        mean_two = 0;
      }
      prob_three = 1 - (prob_one + prob_two);
      if (prob_three > 0)
      {
        mean_three = (cum_mean[255] - cum_mean[i]) / prob_three;
      }else
      {
        mean_three = 0;
      }
      in_btw_var[k][i] = prob_one * (mean_one - global_av) * (mean_one - global_av) +  prob_two * (mean_two - global_av) * (mean_two - global_av) +  
      prob_three * (mean_three - global_av) * (mean_three - global_av);
      if (in_btw_var[k][i] > max)
      {
        max = in_btw_var[k][i];
        threshold_one = k;
        threshold_two = i;
      }
    }
  } 
  printf("T1: %f : T2: %f Max: %f \n",threshold_one, threshold_two,max);
  double **pt_new_image = image_create_image_channel(height,width);
  for (int x = 0; x < height; x++)
  {
    for (int y = 0; y < width;y++)
    {
      if (pt_image[x][y] < threshold_one)
      {
        pt_new_image[x][y] = 0.0;
      }else if (pt_image[x][y] < threshold_two)
      {
        pt_new_image[x][y] = 128.0;
      } else 
      {
        pt_new_image[x][y] = 255.0;
      }     
    }
  }
  return (pt_new_image);
}
/**
 * @brief performs variable thresholding using moving averages
 * @param pt_img pointer to an image channel
 * @param height number of rows in the image channel
 * @param width number of columns in the image channel
 * @param n number of points to use in calculating the moving average
 * @param scaler proportion of the average to be used as a threshold
 * @remark Threshold is calculated as scaler * moving average
*/
#define SNEILMAXGRAYVALUE 255.0F
#define SNEILMINGRAYVALUE 0.0F


double** image_seg_moving_average_thres(double **pt_img, size_t height, size_t width, size_t n,double scaler)
{
  size_t startn,num_points;
  double sum,average,threshold,**pt_seg_img;
  pt_seg_img = image_create_image_channel(height,width);
  for (size_t x = 0; x < height; x++)
  {
    startn = 0;
    num_points = 0;
    sum = 0;
    for (size_t y = 0; y < width;y++)
    {
      sum += pt_img[x][y];
      num_points += 1;
      if (num_points > n)
      {
        sum -= pt_img[x][startn];
        startn += 1;
        num_points = n;
      }
      average = sum / num_points;
      threshold = scaler * average;
      if (pt_img[x][y] > threshold)
      {
        pt_seg_img[x][y] = SNEILMAXGRAYVALUE;
      }else
      {
        pt_seg_img[x][y] = SNEILMINGRAYVALUE;
      }
    }
  }
  return (pt_seg_img);
}
/**
 * @brief performs the k means clustering algorithm
*/
double** image_seg_k_means(double **pt_img, size_t height, size_t width,int k,double threshold)
{
  double *class_mean = (double *)calloc(k, sizeof(double)); // array to hold the class means
  double *class_mean_sum = (double *)calloc(k,sizeof(double)); // array to hold the sum of pixels in for each class
  int *class_mean_count = (int *)calloc(k,sizeof(int)); // array to hold the number of pixels in each class
  int x, y ,istaken = FALSE;
  time_t t;
  printf("Initial Means: ");
  for (int i = 0; i < k; i++)
  {
    do
    {
      srand((unsigned)time(&t));
      istaken = FALSE;
      x = rand() % height;
      y = rand() % width;
      class_mean[i] = pt_img[x][y];
      for (int j = 0; j < i; j++)
      {
        if (class_mean[j] == class_mean[i])
        {
          istaken = TRUE;
          break;
        } 
      }
    } while (istaken);   
    printf("%3d ",(int)class_mean[i]);   
  }
  printf("\n");
  double **label = image_create_image_channel(height,width); // 2D array to hold a class label for a corresponding pixel location
  double mean_dist,pixel_mean_dist_prev,pixel_mean_dist_cur,pixel_value,new_mean ;
  do
  {
    mean_dist = 0;
    //assign pixels to the nearest mean
    for (size_t x = 0; x < height; x++)
    {
      for (size_t y = 0; y < width;y++)
      {
        pixel_value = pt_img[x][y];
        pixel_mean_dist_prev = (pixel_value - class_mean[0]) * (pixel_value - class_mean[0]);
        label[x][y] = 0.0;
        for (size_t i = 0; i < k; i++)
        {
          pixel_mean_dist_cur = (pixel_value - class_mean[i]) * (pixel_value - class_mean[i]);
          if (pixel_mean_dist_cur < pixel_mean_dist_prev)
          {
            label[x][y] = (double)i;
            pixel_mean_dist_prev = pixel_mean_dist_cur; 
          }    
             
        }
        class_mean_sum[(int)label[x][y]] += pixel_value;
        class_mean_count[(int)label[x][y]] += 1;
      }
    }
    //update the mean
    for (size_t i = 0; i < k; i++)
    {
      new_mean = class_mean_sum[i] / (double)class_mean_count[i];     
      mean_dist += (class_mean[i] - new_mean) * (class_mean[i] - new_mean);
      class_mean[i] = new_mean;
      class_mean_sum[i] = 0;
      class_mean_count[i] = 0;
    }  
  } while (mean_dist > threshold);
  printf("Final Means: ");
  double value = 255 / (k - 1);
  for (size_t i = 0; i < k; i++)
  {
    printf("%3d ",(int)class_mean[i]);
  }
  printf("\n");
  for (size_t x = 0; x < height; x++)
  {
    for (size_t y = 0; y < width;y++)
    {
      label[x][y] = round(label[x][y] * value);
    }
  }
  free(class_mean);
  free(class_mean_count);
  free(class_mean_sum);
  return (label);
}

/**
 * @brief implements SLIC algorithm
*/

double** image_seg_slic(double **pt_img,int height, int width, int num_super_pixel, double threshold, double scaler)
{
  int num_tp, sp_interval,num_spy, num_spx,sp_cx,sp_cy;
  num_tp = height * width; // total number of pixels in the image
  sp_interval = sqrt(((double)(num_tp)) / ((double)num_super_pixel)); // grid spacing
  num_spx = height / sp_interval;
  num_spy = width / sp_interval;
  sp_cx = sp_interval / 2;
  sp_cy = sp_cx;
  int p = 0;
  printf("ntp: %d : s: %d : nspx: %d : nspy: %d\n",num_tp,sp_interval,num_spx,num_spy);
  typedef struct means
  {
    double intensity;
    int pos_x;
    int pos_y;
    double intes_sum;
    int pos_x_sum;
    int pos_y_sum;
    double count;

  } MEANS, *PMEANS, **PPMEANS;
  typedef struct superlabel
  {
    int mean_x;
    int mean_y;
    double_t distance;
  } SUPL, *PSUL, **PPSUL;
  PPMEANS pt_sp_means = (PPMEANS)malloc(sizeof(PMEANS) * num_spx);
  for (size_t i = 0; i < num_spx; i++)
  {
    pt_sp_means[i] = calloc(num_spy, sizeof(MEANS));
  }
  int cx, cy, tx,ty;
  double gradientx, gradienty, gradcur, gradprev;
  for (size_t x = 0; x < num_spx; x++)
  {
    cx = x * sp_interval + sp_cx;
    for (size_t y = 0; y < num_spy; y++)
    {
      cy = y * sp_interval + sp_cy;
      pt_sp_means[x][y].pos_x = cx;
      pt_sp_means[x][y].pos_y = cy;
      pt_sp_means[x][y].intensity = pt_img[x][y];
      gradprev = 255.0;
      // calculate the gradient in a 3 x 3 neighborhood about this center
      for (int s = -1; s <= 1; s++)
      {
        tx = cx + s;        
        for (int t = -1; t <= 1; t++)
        {
          ty = cy + t;
          gradientx = pt_img[tx + 1][ty] -  pt_img[tx][ty];
          gradientx = gradientx < 0 ? -gradientx : gradientx;
          gradienty = pt_img[tx][ty + 1] -  pt_img[tx][ty];
          gradienty = gradienty < 0 ? -gradienty : gradienty;
          gradcur = gradientx + gradienty;
          if (gradcur < gradprev)
          {
            gradprev = gradcur;
            pt_sp_means[x][y].pos_x = tx;
            pt_sp_means[x][y].pos_y = ty;
            pt_sp_means[x][y].intensity = pt_img[tx][ty];
          }
        }
      }  
    }   
  }
  printf("Done 1\n");

  PPSUL labels = (PPSUL)malloc(sizeof(PSUL) * height);
  for (size_t i = 0; i < height; i++)
  {
    labels[i] = malloc(sizeof(SUPL) * width);
    if (labels[i] == NULL)
    {
      printf("Failed: \n");
    }    
    for (size_t j = 0; j < width; j++)
    {
      labels[i][j].distance = 99999999999999;
      labels[i][j].mean_x = -1;
      labels[i][j].mean_y = -1;
    }  
  }
  printf("Done 2\n");
  double intes, dc,dx,dy,dist, residual_error;
  int startx, starty, endx, endy,ind_x, ind_y;
  MEANS new_mean;
  scaler = (scaler * scaler) / ((double)(sp_interval * sp_interval));
  do
  {
    for (size_t x = 0; x < num_spx; x++)
    {
      for (size_t y = 0; y < num_spy; y++)
      {
        intes = pt_sp_means[x][y].intensity;
        cx = pt_sp_means[x][y].pos_x;
        cy = pt_sp_means[x][y].pos_y;
        startx = cx - sp_interval;
        starty = cy - sp_interval;
        endx = cx + sp_interval;
        endy = cy + sp_interval;
        if (startx < 0)
        {
          startx = 0;
        }
        if (starty < 0)
        {
          starty = 0;
        }
        if (endx >= height)
        {
          endx = height-1;
        }
        if (endy >= width)
        {
          endy = width - 1;
        }
        for (size_t i = startx; i <= endx; i++)
        {
          dx = (pt_sp_means[x][y].pos_x - i) *  (pt_sp_means[x][y].pos_x - i);
          for (size_t j = starty; j <= endy; j++)
          {
            dy = (pt_sp_means[x][y].pos_y - j) *  (pt_sp_means[x][y].pos_y - j);
            dc =  (pt_img[i][j] - intes) * (pt_img[i][j]  - intes) ;
            dist = sqrt((dx + dy)* scaler + dc);
            if (dist < labels[i][j].distance)
            {
              labels[i][j].distance = dc;
              labels[i][j].mean_x = x;
              labels[i][j].mean_y = y;
            }           
          }       
        }
      }
    }
    for (size_t x = 0; x < height; x++)
    {
      for (size_t y = 0; y < width; y++)
      {
        ind_x = labels[x][y].mean_x;
        ind_y = labels[x][y].mean_y;
        if (ind_x < 0 || ind_y < 0)
        {
          continue;
        }        
        pt_sp_means[ ind_x][ind_y].count += 1;
        pt_sp_means[ ind_x][ind_y].intes_sum += pt_img[x][y];
        pt_sp_means[ ind_x][ind_y].pos_x_sum += x;
        pt_sp_means[ ind_x][ind_y].pos_y_sum += y;
      }
    }
    //update the means
    residual_error = 0;
    for (size_t x = 0; x < num_spx; x++)
    {
      for (size_t y = 0; y < num_spy; y++)
      {
        new_mean.intensity = pt_sp_means[x][y].intes_sum / pt_sp_means[x][y].count;
        new_mean.pos_x = pt_sp_means[x][y].pos_x_sum / pt_sp_means[x][y].count;
        new_mean.pos_y = pt_sp_means[x][y].pos_y_sum / pt_sp_means[x][y].count;
        residual_error += (new_mean.intensity - pt_sp_means[x][y].intensity) * (new_mean.intensity - pt_sp_means[x][y].intensity) + (new_mean.pos_x - pt_sp_means[x][y].pos_x) * (new_mean.pos_x - pt_sp_means[x][y].pos_x) + (new_mean.pos_y - pt_sp_means[x][y].pos_y) * (new_mean.pos_y - pt_sp_means[x][y].pos_y);
        pt_sp_means[x][y].count = 0;
        pt_sp_means[x][y].intensity = new_mean.intensity;
        pt_sp_means[x][y].intes_sum = 0;
        pt_sp_means[x][y].pos_x = new_mean.pos_x;
        pt_sp_means[x][y].pos_x_sum = 0;
        pt_sp_means[x][y].pos_y = new_mean.pos_y;
        pt_sp_means[x][y].pos_y_sum = 0;
      }
    }
    p++;
  } while (residual_error > threshold);
  printf("P = %d\n",p);
 /*for (size_t x = 0; x < num_spx; x++)
  {
    for (size_t y = 0; y < num_spy; y++)
    {
      printf("%6.2f(%3d,%3d) ", pt_sp_means[x][y].intensity, pt_sp_means[x][y].pos_x,pt_sp_means[x][y].pos_y);
    }
    printf("\n");
  }*/
  double **new_img = image_create_image_channel(height,width);
  for (size_t x = 0; x < height; x++)
  {
    for (size_t y = 0; y < width; y++)
    {
      ind_x = labels[x][y].mean_x;
      ind_y = labels[x][y].mean_y;
      if (ind_x < 0 || ind_y < 0)
      {
        new_img[x][y] = 0;
      } else
      {
        new_img[x][y] = pt_sp_means[ind_x][ind_y].intensity;
      }
    }
  }
  for (size_t i = 0; i < num_spx ; i++)
  {
    free(pt_sp_means[i]);
  }
  for (size_t i = 0; i < height; i++)
  {
    free(labels[i]);
  }
  return (new_img);        
}
/**
 * @brief computes the accumulative difference image
*/
void image_seg_adi(double **reference, double **frame,double **adi,int height, int width, double threshold)
{
  double diff;
  for (size_t x = 0; x < height; x++)
  {
    for (size_t y = 0; y < width; y++)
    {
      diff = reference[x][y] - frame[x][y];
      diff = diff < 0 ? -diff : diff;
      adi[x][y] += diff > threshold ? 1 : 0.0;
    }
  }
  return;
}

/**
 * @brief computes the positive accumulative difference image
*/
void image_seg_positive_adi(double **reference, double **frame,double **positive_adi,int height, int width, double threshold)
{
  double diff;
  for (size_t x = 0; x < height; x++)
  {
    for (size_t y = 0; y < width; y++)
    {
      diff = reference[x][y] - frame[x][y];
      positive_adi[x][y] += diff > threshold ? 1 : 0.0;
    }
  }
  return;
}

/**
 * @brief computes the positive accumulative difference image
*/
void image_seg_negative_adi(double **reference, double **frame,double** negative_adi,int height, int width, double threshold)
{
  double diff;
  double neg_threshold = -threshold;
  for (size_t x = 0; x < height; x++)
  {
    for (size_t y = 0; y < width; y++)
    {
      diff = reference[x][y] - frame[x][y];
      negative_adi[x][y] += diff < neg_threshold ? 1 : 0.0;
    }
  }
  return;
}