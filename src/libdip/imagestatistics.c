
double **image_get_unnormalized_histogram(double ***image, size_t image_height, size_t image_width, int num_of_channels)
{
  double **stat = (double **)malloc(sizeof(double *) * num_of_channels);
  for (size_t n = 0; n < num_of_channels; n++)
  {
    stat[n] = (double *)calloc(256,sizeof(double));
    image_sum_intensity(image[n],stat[n],image_height,image_width);

  }
  return (stat);

}
double **image_get_normalized_histogram(double ***image, size_t image_height, size_t image_width, int num_of_channels)
{
  double **stat = (double **)malloc(sizeof(double *) * num_of_channels);
  double divisor = image_height * image_width;
  int r;
  for (size_t n = 0; n < num_of_channels; n++)
  {
    stat[n] = (double *)calloc(256,sizeof(double));

      for (size_t x = 0; x < image_height; x++)
      {
        for (size_t y = 0; y < image_width; y++)
        {
          r = (int)image[n][x][y];
          stat[n][r] += 1.0F;

        }
    
      }
      for (size_t i = 0; i < 256; i++)
      {
        stat[n][i] = stat[n][i] / divisor;
      }
      
  }
  return (stat);
}
double **image_get_localnormalized_histogram(double ***image, size_t orig_x, size_t orig_y,size_t end_x,size_t end_y,double dim,int num_of_channels)
{
  double **stat = (double **)malloc(sizeof(double *) * num_of_channels);
  int r;
  for (size_t n = 0; n < num_of_channels; n++)
  {
    stat[n] = (double *)calloc(256,sizeof(double));

    for (size_t x = orig_x; x <= end_x; x++)
    {
      for (size_t y = orig_y; y <= end_y; y++)
      {
        r = (int)image[n][x][y];
        stat[n][r] += 1.0F;

      }
  
    }
    for (size_t i = 0; i < 256; i++)
    {
      stat[n][i] = stat[n][i] / dim;
    }
      
  }
  return (stat);
}
void image_sum_intensity(double **image_channel,double *arr ,size_t image_height, size_t image_width)
{
  int r;
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      r = (int)image_channel[x][y];
      arr[r] += 1.0F;

    }
    
  }
  
  
  return;
}

void image_get_mean(double ***image,double **stat,double *mean_var,int num_of_channels)
{
  double mean;
  for (size_t n = 0; n < num_of_channels; n++)
  {
    mean = 0.0F;
    for (size_t i = 0; i < 256; i++)
    {
      mean += stat[n][i] * ((double)(i));
    }
    *(mean_var + n * 2 ) = mean;
  }
}

void image_get_sample_mean(double ***image,size_t image_height,size_t image_width,int num_of_channels,double *mean_var)
{
  double sample_mean;
  double sample_space = image_height * image_width;
  for (size_t n = 0; n < num_of_channels; n++)
  {
    sample_mean = 0.0F;
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        sample_mean += image[n][x][y];
      }
      
    }
    *(mean_var + n * 2 ) = sample_mean / sample_space;
  }
  return;
}


double *image_get_nmoment(double ***image,double **stat,int num_of_channels, int p)
{
  double diff;
  double var;
  double *mean_var = (double *)malloc(sizeof(double) * (num_of_channels * 2));
  image_get_mean(image,stat,mean_var,num_of_channels);
  for (size_t n = 0; n < num_of_channels; n++)
  {
    var = 0.0F;
    for (size_t i = 0; i < 256; i++)
    {
      diff = i - (*(mean_var + n * 2));
      var += pow(diff,p) * stat[n][i];

    }
    *(mean_var + n * 2 + 1) = sqrt(var);
  } 
  return (mean_var);
}
double *image_get_local_mean_moment(double ***image,size_t image_height, size_t image_width, int num_of_channels, int p, size_t cx, size_t cy, int neig_x, int neig_y)
{
  int divisorx = neig_x /2;
  int divisory = neig_y / 2;
  int orig_x = cx - divisorx;
  int orig_y = cy - divisory;                            
  int end_x = cx + divisorx;
  int end_y = cy + divisory;
  
  orig_x = orig_x < 0? 0: orig_x;
  orig_y = orig_y < 0? 0: orig_y;
  end_x = end_x >= image_height ? image_height - 1 : end_x;
  end_y = end_y >= image_width ? image_width - 1 : end_y;
  double dim = neig_x * neig_y;
  double **stat = image_get_localnormalized_histogram(image,orig_x,orig_y,end_x,end_y,dim,num_of_channels);
  double *mean_var = image_get_nmoment(image,stat,num_of_channels,p);
  return (mean_var);
}

double ***image_local_histogram_en(double ***image,size_t image_height, size_t image_width,int num_of_channels,int neig_x,int neig_y,double c, double k_0, double k_1,double k_2, double k_3)
{
  double ***new_image = image_create_image(image_height,image_width,num_of_channels);
  
  double **stat= image_get_normalized_histogram(image,image_height,image_width,num_of_channels);
  double *global_mean_var = image_get_nmoment(image,stat,num_of_channels,2);
  printf("Global mean: %f\n",(*(global_mean_var + 0)));
  printf("Local Variance: %f\n",(*(global_mean_var + 1)));
  double *local_mean_var;
  int candidate;
  double global_mean;
  double local_mean;
  double global_var;
  double local_var;
  for (size_t n = 0; n < num_of_channels; n++)
  {
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        local_mean_var = image_get_local_mean_moment(image,image_height,image_width,num_of_channels,2,x,y,neig_x,neig_y);
        global_mean = (*(global_mean_var + n * 2));
        global_var = *(global_mean_var + n * 2 + 1);
        local_mean = *(local_mean_var + n * 2);
        local_var = *(local_mean_var + n * 2 + 1);
        candidate = (((k_0 * global_mean) <= local_mean ) && (local_mean <= (k_1 * global_mean))) && (((k_2 * global_var) <= local_var ) && (local_var <= (k_3 * global_var)));
        if (candidate)
        {
          new_image[n][x][y] = c * image[n][x][y];
        } else
        {
          new_image[n][x][y] = image[n][x][y];
        }
        free(local_mean_var);
      }
      
    }
    
  }
  return (new_image);
  
}

double** image_local_histogram_equalization(double **image,size_t height, size_t width, int m, int n)
{
  double** img = image_create_image_channel(height,width);
  double val[256];
  double sum = 0;
  for (size_t i = 0; i < 256; i++)
  {
    val[i] = 0.0;
  }
  double scaler = 1.0  /((double)(m * n));
  int a = m / 2, b = n / 2, tx , ty,xbound = height + height - 2, ybound = width + width - 2;
  for (size_t x = 0; x < height; x++)
  {
    for (size_t y = 0; y < width; y++)
    {
      for (int s = -a; s <= a; s++)
      {
        tx = x + s;
        if (tx < 0)
        {
        tx = x - s;
        }else if (tx >= height)
        {
            tx = xbound - tx;
        }       
        for (int t = -b; t <= b; t++)
        {
          ty = y + t;
          if (ty < 0)
          {
              ty = y - t;
          }else if (ty >= width)
          {
                ty = ybound - ty;
          } 
          val[(int)image[tx][ty]] += scaler;
        }
        
      }
      for (size_t i = 0; i <= ((int)(image[x][y])); i++)
        {
          sum += val[i];
        }
        img[x][y] = 255.0 * sum;
        sum = 0.0;
        for (size_t i = 0; i < 256; i++)
        {
          val[i] = 0.0;
        }    
    } 
  }
  return (img);
}