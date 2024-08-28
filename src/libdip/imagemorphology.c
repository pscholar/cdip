//#define CAST2D(l) ((int (*)[])l)
#define CAST2D(l) ((double (*)[])l)
#define IDONTCARE 600
typedef int BOOLEAN;
int counter;
/**
 * @brief performs grayscale erosion of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @returns eroded image
*/
double **image_morphology_grayscale_erosion(double **image,size_t image_height, size_t image_width, int rowsB, int colsB)
{
  double **eroded_image = image_create_image_channel(image_height,image_width);
  int tx,ty,a,b;
  int isNotValid;
  double min;
  a = rowsB / 2;
  b = colsB / 2;
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        min = image[x][y];
        for (int s = -a; s <=  a; s++)
        {
          tx = x + s;
          for (int t = 0; t <= b; t++)
          {
            ty = y + t;
            isNotValid = (tx >= image_height) || (tx < 0) || (ty >= image_width) || (ty < 0);
            if (!isNotValid)
            {
              min = image[tx][ty] < min ? image[x][y] : min;
            }            
          }          
        }
        eroded_image[x][y] = min;
      }      
    }    
  return (eroded_image);
}

/**
 * @brief performs grayscale dilation of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @returns dilated image
*/
double **image_morphology_grayscale_dilation(double **image,size_t image_height, size_t image_width, int rowsB, int colsB)
{
  double **eroded_image = image_create_image_channel(image_height,image_width);
  int tx,ty,a,b;
  int not_valid;
  double max;
  a = rowsB / 2;
  b = colsB / 2;
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      max = image[x][y];
      for (int s = -a; s <=  a; s++)
      {
        tx = x - s;
        for (int t = 0; t <= b; t++)
        {
          ty = y - t;
          not_valid = (tx >= image_height) || (tx < 0) || (ty >= image_width) || (ty < 0);
          if (not_valid)
          {
            max = image[tx][ty] > max ? image[x][y] : max;
          }            
        }          
      }
      eroded_image[x][y] = max;
    }      
  }    
  return (eroded_image);
}

/**
 * @brief performs grayscale opening of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @returns opened image
*/
double **image_morphology_grayscale_opening(double **image,size_t image_height, size_t image_width,int rowsB, int colsB)
{
  double **eroded_image = image_morphology_grayscale_erosion(image,image_height,image_width,rowsB,colsB);
  double **opening = image_morphology_grayscale_dilation(eroded_image,image_height,image_width,rowsB,colsB);
  image_free_kernel(eroded_image,image_height);
  return (opening);
}

/**
 * @brief performs grayscale closing of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @returns closed image
*/
double **image_morphology_grayscale_closing(double **image,size_t image_height, size_t image_width,int rowsB, int colsB)
{
  double **dilated_image = image_morphology_grayscale_dilation(image,image_height,image_width,rowsB,colsB);
  double **closed= image_morphology_grayscale_erosion(dilated_image,image_height,image_width,rowsB,colsB);
  image_free_kernel(dilated_image,image_height);
  return (closed);
}

/**
 * @brief performs grayscale morphological smoothing of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @returns smoothed image
*/
double **image_morphology_grayscale_smoothing(double **image,size_t image_height, size_t image_width,int rowsB, int colsB)
{
  double **opened_image = image_morphology_grayscale_opening(image,image_height,image_width,rowsB,colsB);
  double **smoothed_image = image_morphology_grayscale_closing(opened_image,image_height,image_width,rowsB,colsB);
  image_free_kernel(opened_image,image_height);
  return (smoothed_image);
}

/**
 * @brief performs grayscale morphological gradient of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @returns morphology gradient image
*/
double **image_morphology_grayscale_gradient(double **image,size_t image_height, size_t image_width, int rowsB, int colsB)
{
  double **dilated_image = image_morphology_grayscale_dilation(image,image_height,image_width,rowsB,colsB);
  double **eroded_image = image_morphology_grayscale_erosion(image,image_height,image_width,rowsB,colsB);
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      dilated_image[x][y] -= eroded_image[x][y];
    }      
  }    
  
  image_free_kernel(eroded_image,image_height);
  return (dilated_image);
}

/**
 * @brief performs grayscale morphological top hat transform of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @returns top hat transform of image
*/
double **image_morphology_grayscale_tophat_transform(double **image,size_t image_height, size_t image_width, int rowsB, int colsB)
{
  double **opened_image = image_morphology_grayscale_opening(image,image_height,image_width,rowsB,colsB);

  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      opened_image[x][y] = image[x][y] - opened_image[x][y];
    }      
  }    
  return (opened_image);
}

/**
 * @brief performs grayscale morphological bottom hat transform of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @returns bottom hat transform of image
*/
double **image_morphology_grayscale_bottomhat_transform(double **image,size_t image_height, size_t image_width,int rowsB, int colsB)
{
  double **closed_image = image_morphology_grayscale_closing(image,image_height,image_width,rowsB,colsB);
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      closed_image[x][y] -= image[x][y];
    }      
  }    
  return (closed_image);
}

/**
 * @brief performs grayscale morphological granulometry of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param m maximum size of structuring element to be used(odd number)
 * @returns pointer to an arrays containing differences between openings 
*/
double * image_morphology_grayscale_granulometry(double **image,size_t image_height, size_t image_width,int m)
{
  double **opened_image;
  double *ptr = (double *)calloc(m + 1, sizeof(double));
  double sum;
  int size;
  int end = m / 2;
  for (int i = 0; i <= end; i++)
  {
    size = 2 * i + 1;
    opened_image = image_morphology_grayscale_opening(image,image_height,image_width,size,size);
    sum = 0;
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        sum += opened_image[x][y];
      }      
    }
    ptr[size - 1] = sum;     
    image_free_kernel(opened_image,image_height);    
  } 
  return (ptr);
}

/**
 * @brief performs grayscale morphological geodesic dilation of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @param m size of the geodesic dilaition
 * @returns geodesic dilation of the image
*/
double **image_morphology_grayscale_geodesic_dilation(double **image,size_t image_height, size_t image_width,int rowsB, int colsB,int m)
{
  int i = 0;
  double **dilation_new;
  double **dilation_old;
  if(m < 1) return (NULL);
  dilation_old = image_morphology_grayscale_dilation(image,image_height,image_width,rowsB,colsB);
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      dilation_old[x][y] = dilation_old[x][y] > image[x][y] ? image[x][y]:dilation_old[x][y];
    }      
  } 
  i++;
  while (i < m)
  {
    dilation_new = image_morphology_grayscale_dilation(dilation_old,image_height,image_width,rowsB,colsB);
    image_free_kernel(dilation_old,image_height);
    dilation_old = dilation_new;
    for (size_t x = 0; x < image_height; x++)
    {
      for (size_t y = 0; y < image_width; y++)
      {
        dilation_old[x][y] = dilation_old[x][y] > image[x][y] ? image[x][y]:dilation_old[x][y];
      }      
    }    
    i++;
  }
  return (dilation_old);
}


/**
 * @brief performs grayscale morphological reconstruction by geodesic dilation of an image
 * @param image pointer to image
 * @param image_height number of rows in image
 * @param image_width number of columns in the image
 * @param num_of_channels number of channels in the image
 * @param rowsB number of rows in structuring element
 * @param colsB number of columns in structuring element
 * @param m size of the geodesic dilaition
 * @returns reconstructed of the image
*/
double **image_morphology_grayscale_dilation_reconstruction(double **image,size_t image_height, size_t image_width,int rowsB, int colsB,int m)
{
  double **dilated_old = image_morphology_grayscale_geodesic_dilation(image,image_height,image_width,rowsB,colsB,1);
  int isSame = image_compare(image,dilated_old,image_height,image_width);
  double **dilated_new;
  while (!isSame)
  {
    dilated_new = image_morphology_grayscale_geodesic_dilation(dilated_old,image_height,image_width,rowsB,colsB,1);
    isSame = image_compare(dilated_new,dilated_old,image_height,image_width);
    image_free_kernel(dilated_old,image_height);
    dilated_old = dilated_new;  
  }  
  return (dilated_old);
}

/**
 * @brief compares two images
*/
int image_compare(double **image1,double **image2,size_t image_height, size_t image_width)
{
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      if(((int)image1[x][y]) != ((int)image2[x][y]))
      {
        return (0);
      }
    }      
  } 
  return (1);
}

/**
 * @brief finds the subtraction of a set
 */
void image_subtraction(int image_height, int image_width,double **image1,double **image2)
{
  
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      image2[x][y] = image1[x][y] - image2[x][y]; 
    }      
  } 

  return;
}

/**
 * @brief finds the intersection of a set
 */
void  image_intersection(int image_height, int image_width,double **image1,double **image2)
{
  
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      if((int)image1[x][y] != (int)image2[x][y])
      {
        image2[x][y] = 0.0F;
      } 
    }      
  } 
  return;
}

/**
 * @brief copies elements of an array into another
 */
void image_copy_elements(int image_height, int image_width,double **image1,double **image2)
{
  
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      image2[x][y] = image1[x][y];
    }      
  } 
  
}

/**
 * @brief unions
 */
void image_union(int image_height, int image_width,double **image1,double **image2)
{
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      if ((int)image1[x][y] > 0)
      {
        image2[x][y] = image1[x][y];
      }    
    }      
  } 
}

/**
 * @brief checks whether an image consists of only zeros.
*/
int image_is_zeroed(int image_height, int image_width,double **image)
{
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      if ((int)image[x][y])
      {
        return FALSE;
      }    
    }      
  }  
  return TRUE;
}
/**
 * @brief get the complement of an image
*/
void image_complement(int image_height, int image_width,double **image)
{
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      image[x][y] = 255.0 - image[x][y];
    }      
  }  
  return;
}
/**
 * @brief performs  binary erosion
*/

double **image_binary_erosion(double **image,size_t image_height, size_t image_width, int m , int n, double (*B)[n])
{
  double **new_image = image_create_image_channel(image_height,image_width);
  int a,b,tx,ty, valB, valA;
  BOOLEAN isSubSet = TRUE, isOutofRange;
  a = m / 2; b = n / 2;
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width; y++)
    {
      for (int s = -a; s <= a; s++)
      {
        if (!isSubSet)
        {
          break;
        }  
        tx = x + s;              
        for (int t = -b; t <= b; t++)
        {
          ty = y + t;
          valB = (int)B[s + a][t + b];
          isOutofRange = (tx < 0) || (ty < 0) || (tx >= image_height) || (ty >= image_width);
          if (isOutofRange)
          {
            /*Assume zero padding of the array*/           
            if( valB != IDONTCARE)
            {
              isSubSet = FALSE;
              break;
            }
          }else
          {
            valA = (int)image[tx][ty];
            if ((valB != IDONTCARE) && (valA !=255))
            {
              isSubSet = FALSE;
              break;
            }            
          }                   
        }        
      }
      if (isSubSet)
      {
        new_image[x][y] = 255.0F;
      }
      isSubSet = TRUE;
    }    
  }  
  return (new_image);
}

/**
 * @brief performs binary dilation
*/
double **image_binary_dilation(double **image,size_t image_height, size_t image_width, int m , int n, double (*B)[n])
{
  double **new_image = image_create_image_channel(image_height,image_width);
  int a,b,tx,ty, valB, valA;
  BOOLEAN isSubSet = FALSE, isOutofRange;
  a = m / 2; b = n / 2;
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width; y++)
    {
      for (int s = -a; s <= a; s++)
      { 
        if (isSubSet)
        {
          break;
        }  
        tx = x + s;              
        for (int t = -b; t <= b; t++)
        {
          ty = y + t;
          valB = (int)B[s + a][t + b];
          isOutofRange = (tx < 0) || (ty < 0) || (tx >= image_height) || (ty >= image_width);
          if (!isOutofRange)
          {
            valA = (int)image[tx][ty];
            if ((valA == 255) & (valB == 255))
            {
              isSubSet = TRUE;
              break;
            }      
          }               
        }        
      }
      if (isSubSet)
      {
        new_image[x][y] = 255.0F;
      }
      isSubSet = FALSE;
    }    
  }
  
  return (new_image);
}

/**
 * @brief performs binary opening
*/
double **image_binary_opening(double **image,size_t image_height, size_t image_width, int m , int n, double (*B)[n])
{
  double **new_image = image_create_image_channel(image_height,image_width);
  int a,b,*ptr,tx,ty, valB, valA;
  BOOLEAN isSubSet = TRUE, isOutofRange;
  a = m/ 2; b = n / 2;

  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width; y++)
    {
      for (int s = -a; s <= a; s++)
      {
        if (!isSubSet)
        {
          break;
        }  
        tx = x + s;              
        for (int t = -b; t <= b; t++)
        {
          ty = y + t;
          valB = (int)B[s + a][t + b];
          isOutofRange = (tx < 0) || (ty < 0) || (tx >= image_height) || (ty >= image_width);
          if (isOutofRange)
          {
            /*Assume zero padding of the array*/           
            if( valB != IDONTCARE)
            {
              isSubSet = FALSE;
              break;
            }
          }else
          {
            valA = (int)image[tx][ty];
            if ((valB != IDONTCARE) && (valA != 255))
            {
              isSubSet = FALSE;
              break;
            }            
          }                   
        }        
      }
      if (isSubSet)
      {
        for (int s = -a; s <= a; s++)
        { 
          tx = x + s;              
          for (int t = -b; t <= b; t++)
          {
            ty = y + t;
            isOutofRange = (tx < 0) || (ty < 0) || (tx >= image_height) || (ty >= image_width);
            if (!isOutofRange)
            {
              new_image[ty][ty] = 255.0F;              
            }                
          }        
        }       
      }
      isSubSet = TRUE;
    }    
  }      
  
  
  
  return (new_image);
}

/**
 * @brief performs binary closing
*/
double **image_binary_closing(double **image,size_t image_height, size_t image_width, int m , int n, double (*B)[n])
{
   double **new_image = image_create_image_channel(image_height,image_width);
  int a,b,*ptr,tx,ty, valB, valA;
  BOOLEAN isSubSet = FALSE, isOutofRange;
  a = m/ 2; b = n/ 2;
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width; y++)
    {
      for (int s = -a; s <= a; s++)
      {
        if (isSubSet)
        {
          break;
        }  
        tx = x + s;              
        for (int t = -b; t <= b; t++)
        {
          ty = y + t;
          valB = B[s + a][t + b];
          isOutofRange = (tx < 0) || (ty < 0) || (tx >= image_height) || (ty >= image_width);
          if (!isOutofRange)
          {
            valA = (int)image[tx][ty];
            if ((valA == 255) & (valB == 255))
            {
              isSubSet = TRUE;
              break;
            }     
          }               
        }        
      }
      if (!isSubSet)
      {
        for (int s = -a; s <= a; s++)
        { 
          tx = x + s;              
          for (int t = -b; t <= b; t++)
          {
            ty = y + t;
            isOutofRange = (tx < 0) || (ty < 0) || (tx >= image_height) || (ty >= image_width);
            if (!isOutofRange)
            {
              new_image[x][y] = 255.0F;              
            }                
          }        
        }  
      }
      isSubSet = FALSE;
    }    
  }      /* code */

  image_complement(image_height,image_width,new_image);
  return (new_image);
}

/**
 * @brief performs a binary hit or mis transform
*/
double **image_binary_hit_or_miss(double **image,size_t image_height, size_t image_width, int m , int n, double (*B)[n])
{
  double **new_image = image_create_image_channel(image_height,image_width);
  int a,b,tx,ty, valB, valA;
  BOOLEAN isSubSet = TRUE, isOutofRange;
  a = m / 2; b = n / 2;
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width; y++)
    {
      for (int s = -a; s <= a; s++)
      {
        if (!isSubSet)
        {
          break;
        }  
        tx = x + s;              
        for (int t = -b; t <= b; t++)
        {
          ty = y + t;
          valB = (int)B[s + a][t + b];
          isOutofRange = (tx < 0) || (ty < 0) || (tx >= image_height) || (ty >= image_width);
          if (isOutofRange)
          {
            /*Assume zero padding of the array*/           
            if( (valB != IDONTCARE) && ( valB != 0))
            {
              isSubSet = FALSE;
              break;
            }
          }else
          {
            valA = (int)image[tx][ty];
            if ((valB != IDONTCARE) && (valA != valB))
            {
              isSubSet = FALSE;
              break;
            }            
          }                   
        }        
      }
      if (isSubSet)
      {
        new_image[x][y] = 255.0F;
      }
      isSubSet = TRUE;
    }    
  } 
  return (new_image);
}

/**
 * @brief performs binary morphological boundary extraction
*/
double **image_boundary_morphology(double **image,size_t image_height, size_t image_width, int m , int n, double (*B)[n])
{
  double **new_image = image_binary_erosion(image,image_height,image_width,m,n,B);
  for (size_t x = 0; x < image_height; x++)
  {
    for (size_t y = 0; y < image_width; y++)
    {
      new_image[x][y] = image[x][y] - new_image[x][y];
    }
  }  
  return (new_image);
}

/**
 * @brief extracts connected
*/
double **image_binary_extract_connected_components(double **image,size_t image_height, size_t image_width,int m , int n, double (*B)[n], int connPosX, int connPosY)
{
  double **new_image = image_create_image_channel(image_height,image_width);
  new_image[connPosX][connPosY] = 255.0F;
  double **p_image =  image_binary_dilation(new_image,image_height, image_width,m,n,B);
  image_intersection(image_height,image_width,image,new_image);
  int is_same = image_compare(new_image,p_image,image_height,image_width);
  image_free_kernel(new_image,image_height);
  new_image = p_image;
  while (!is_same)
  {
    p_image =  image_binary_dilation(new_image,image_height, image_width,m,n,B);
    image_intersection(image_height,image_width,image,new_image);
    is_same = image_compare(new_image,p_image,image_height,image_width);
    image_free_kernel(new_image,image_height);
    new_image = p_image;
  }
  image_union(image_height,image_width,image,new_image);
  return (new_image);
}

/**
 * @brief performs the morphological thining.
*/

double** image_binary_thining(double **image,size_t image_height, size_t image_width,int m , int n,int count, double *B[count])
{
  BOOLEAN noChangePass = TRUE;
  BOOLEAN isCheckDone = FALSE;
  double **prev_image = image_create_image_channel(image_height,image_width);
  image_copy_elements(image_height,image_width,image,prev_image);
  double **cur_image;
  counter = 0;
  do
  {    
    for (int i = 0; i < count; i++)
    {
      cur_image = image_binary_hit_or_miss(prev_image,image_height, image_width, m,n,CAST2D(B[i]));
      image_subtraction(image_height,image_width,prev_image,cur_image);
      if (!isCheckDone)
      {
        noChangePass = image_compare(prev_image,cur_image,image_height,image_width);
        if (!noChangePass)
        {
          isCheckDone = TRUE;
        } 
      }
      image_free_kernel(prev_image,image_height);
      prev_image = cur_image;
    }
    counter++;
    isCheckDone = FALSE;
  } while (!noChangePass);
  return (cur_image); 
}

/**
 * @brief performs the binary skeleton.
*/

double** MorphologicalBinarySkeleton(double **image,size_t image_height, size_t image_width,int m , int n, double (*B)[n])
{
  counter = 0;
  double **skeleton = image_create_image_channel(image_height,image_width);
  double **erosion = image_create_image_channel(image_height,image_width);
  image_copy_elements(image_height,image_width,image,erosion);
  double **opening;
  double **cur_erosion;
  while (!image_is_zeroed(image_height,image_width,erosion))
  {
    //perform opening
    opening = image_binary_opening(erosion,image_height,image_width,m,n,B);
    image_subtraction(image_height,image_width,erosion,opening);
    image_union(image_height,image_width,opening,skeleton);
    image_free_kernel(opening,image_height);
    cur_erosion = image_binary_erosion(erosion,image_height,image_width,m,n,B);
    image_free_kernel(erosion,image_height);
    erosion = cur_erosion;
    counter++;
  }
  image_free_kernel(erosion,image_height);
  return (skeleton);
}

/**
 * @brief performs the binary geodesic dilation.
 */
double** image_binary_geodesic_dilation(double **mask,size_t image_height, size_t image_width,double **marker,int m, int n, double (*B)[n],int size)
{
  double **dg_old = image_create_image_channel(image_height,image_width);
  image_copy_elements(image_height,image_width,marker,dg_old);
  image_intersection(image_height,image_width,mask,dg_old);
  int i = 0;
  double **dg_new;
  while (i < n)
  {
    dg_new = image_binary_dilation(dg_old,image_height,image_width,m,n,B);
    image_intersection(image_height,image_width,mask,dg_new);
    image_free_kernel(dg_old,image_height);
    dg_old = dg_new;
    i++;
  }  
  return (dg_old);
}

/**
 * @brief performs binary reconstruction by dilation
*/
double** image_binary_geodesic_dilation_reconstruction(double **mask,size_t image_height, size_t image_width,double **marker,int m, int n, double (*B)[n])
{
  counter =0;
  double **dg_old = image_create_image_channel(image_height,image_width);
  image_copy_elements(image_height,image_width,marker,dg_old);
  image_intersection(image_height,image_width,mask,dg_old);
  double **dg_new;
  int is_same;
  do
  {
    dg_new = image_binary_dilation(dg_old,image_height,image_width,m,n,B);
    image_intersection(image_height,image_width,mask,dg_new);
    image_compare(dg_new,dg_old,image_height,image_width);
    image_free_kernel(dg_old,image_height);
    dg_old = dg_new;
  } while (!is_same); 
  return (dg_old);
}

/**
 * @brief does automatic binary hole filling by reconstruction
*/
double** image_hole_filling_by_reconstruction(double **image,size_t image_height, size_t image_width,int m , int n, double (*B)[n])
{
  double **marker = image_create_image_channel(image_height,image_width);
  double **mask= image_create_image_channel(image_height,image_width);
  image_copy_elements(image_height,image_width,image,mask);
  image_complement(image_height,image_width,mask);
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width; y++)
    {
      if ((x == 0) || (x == image_height- 1) || (y == 0) || (y == image_width - 1))
      {
        marker[x][y] = 255 - image[x][y];     
        
      }      
    }
  }
  double** recon = image_binary_geodesic_dilation_reconstruction(mask,image_height,image_width,marker,m, n,B);
  image_complement(image_height,image_width,recon);
  image_free_kernel(mask,image_width);
  image_free_kernel(marker,image_width);
  return (recon);
}

/**
 * @brief does border clearing by reconstruction
*/
double** MorphologicalBorderClearingByReconstruction(double **image,size_t image_height, size_t image_width,int m , int n, double (*B)[n])
{
  double **marker = image_create_image_channel(image_height,image_width);
  for (int x = 0; x < image_height; x++)
  {
    for (int y = 0; y < image_width; y++)
    {
      if ((x == 0) || (x == image_height - 1) || (y == 0) || (y == image_width - 1))
      {
        marker[x][y] = image[x][y];    
        
      }      
    }
  }
  double** recon = image_binary_geodesic_dilation_reconstruction(image,image_height,image_width,marker,m, n,B);
  image_free_kernel(marker,image_height);
  image_subtraction(image_height, image_width,image,recon);
  return (recon);
}

/**
 * @brief does binary opening by reconstruction
*/
double** image_binary_opening_by_reconstruction(double **image,size_t image_height, size_t image_width,int m , int n, double (*B)[n],int k)
{
  counter =0;
  double **eroded_old = image_create_image_channel(image_height,image_width);
  image_copy_elements(image_height,image_width,image,eroded_old);
  double **eroded_new;
  int i = 0;
  while (i < k)
  {
    eroded_new = image_binary_erosion(eroded_old,image_height,image_width,m,n,B);
    image_free_kernel(eroded_old,image_height);
    eroded_old= eroded_new;
    i++;
  }
  double B2[3][3] =  {{1,1,1},{1,1,1},{1,1,1}};
  double** opening = image_binary_geodesic_dilation_reconstruction(image,image_height,image_width,eroded_old,m, n,B2);
  image_free_kernel(eroded_old,image_height);
  return (opening);
}




