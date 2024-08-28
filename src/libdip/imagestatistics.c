
#include <stdio.h>
#include <stdlib.h>

#define MAX_INTENSITY_LEVEL 255
#define MAX_NCHANNELS 3

double *
img_get_unnormalized_histogram(double *image, int width, int height, 
                                  int nchannels)
{ 
  // TODO : Add guard against NULL pointers and failed memory allocation
  size_t row;
  size_t pix_index;
  size_t row_stride = sizeof(double) * height * width * nchannels;
  size_t byte_size = sizeof(double) * nchannels * (MAX_INTENSITY_LEVEL + 1);
  double *hist = (double *)calloc(byte_size, sizeof(double));
  unsigned char r;
  for (int i = 0; i < height; i++)
  {
    row = i * row;
    for (int j = 0; j < width; j++)
    {
      for (int n = 0; n < nchannels; n++)
      {
        pix_index = row + j * nchannels + n;
        r = (unsigned char)image[pix_index];
        hist[r * nchannels + n] += 1.0F;
      }    
    }  
  }
  return hist;
}

double *img_get_normalized_histogram(double *image, int width, int height, 
                                      int nchannels)
{
  // TODO: guard against NULL pointers and division by zero error.
  double *hist = 
                  img_get_unnormalized_histogram(image,width, height,nchannels);

  double num_samples = (double)(height * width);
  int r;
  for (size_t n = 0; n < nchannels; n++)
  {
    for (size_t i = 0; i < MAX_INTENSITY_LEVEL + 1; i++)
    {
      hist[ i * nchannels + n] /= num_samples;
    }
  }    
  return (hist);
}



/**
 * sums up intensity for each channel in an image
*/

void img_sum_intensity(double *image,double *sum,int width, 
                          int height, int nchannels)
{
  size_t row;
  size_t pix_index;
  size_t row_stride = sizeof(double) * height * width * nchannels;
  for (int i = 0; i < height; i++)
  {
    row = i * row;
    for (int j = 0; j < width; j++)
    {
      for (int n = 0; n < nchannels; n++)
      {
        pix_index = row + j * nchannels + n;
        sum[n] += image[pix_index];
      }    
    }  
  }
  return;
}

void img_get_mean_intensity(double **image, double *mean,int width, 
                                int height, int nchannels)
{
  // TODO: guard against NULL pointers and division by zero error
  double nsamples = width * height;
  img_sum_intensity(image,mean,width,height,nchannels);
  
  for (size_t n = 0; n < nchannels; n++)
  {
    mean[n] /= nsamples;
  }

}



void image_get_nmoment(double ***image,double *moment,int width, int height, 
                        int nchannels, int p)
{
  double mean_deviation;
  double pmean_deviation;
  double mean[MAX_NCHANNELS];
  // TODO: guard against NULL pointers
  double *norm_hist = img_get_normalized_histogram(image,width,height,nchannels);
  img_get_mean_intensity(image, mean,width,height,nchannels);

  for (size_t n = 0; n < nchannels; n++)
  {
    pmean_deviation = 0.0F;
    for (size_t i = 0; i < MAX_INTENSITY_LEVEL + 1; i++)
    {
      mean_deviation = i - mean[n];
      pmean_deviation += pow(mean_deviation,p) * norm_hist[i * nchannels + n];

    }
    moment[n] = sqrt(pmean_deviation);
  } 
  free(norm_hist);
}

typedef struct
{
  int x;
  int y;
}Point2D;

typedef struct 
{
  Point2D top_left;
  Point2D botm_right;
} Rect;

void 
img_get_local_mean(double *image,double *mean,int width,int height,
                    int nchannels)
{

}

double ***image_local_histogram_en(double ***image,size_t image_height, 
size_t image_width,int num_of_channels,int neig_x,int neig_y,double c, 
double k_0, double k_1,double k_2, double k_3)
{
  double ***new_image = image_create_image(image_height,image_width,num_of_channels);
  
  double **stat= image_get_normalized_histogram(image,image_height,
  image_width,num_of_channels);
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
        local_mean_var = image_get_local_mean_moment(image,image_height,
        image_width,num_of_channels,2,x,y,neig_x,neig_y);
        global_mean = (*(global_mean_var + n * 2));
        global_var = *(global_mean_var + n * 2 + 1);
        local_mean = *(local_mean_var + n * 2);
        local_var = *(local_mean_var + n * 2 + 1);
        candidate = (((k_0 * global_mean) <= local_mean ) && (
          local_mean <= (k_1 * global_mean))) && (
            ((k_2 * global_var) <= local_var ) && 
            (local_var <= (k_3 * global_var)));
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

double** image_local_histogram_equalization(double **image,size_t height, 
size_t width, int m, int n)
{
  double** img = image_create_image_channel(height,width);
  double val[256];
  double sum = 0;
  for (size_t i = 0; i < 256; i++)
  {
    val[i] = 0.0;
  }
  double scaler = 1.0  /((double)(m * n));
  int a = m / 2, b = n / 2, tx , ty,xbound = height + height - 2, 
  ybound = width + width - 2;
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