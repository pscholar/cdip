short** image_extract_label_regions(double **pt_image, int height,int width, short *set, int *num_of_objects)
{
    int label = 1, A, B, C, D;
    int Bx, By, Dx, Dy, Cx, Cy;
    short equiv_table[2000] = {0}; // 2000 connected components
    short **label_image = (short **)malloc(sizeof(short *) * height);
    for (size_t i = 0; i < height; i++)
    {
        label_image[i] = (short *)calloc(width, sizeof(short));
    }
    for (int x = 0; x < height; x++)
    {
        Bx = Dx = x - 1;
        Cx = x;       
        for (int y = 0; y < width; y++)
        {
            Dy = Cy = y - 1;
            By = y;
            A = pt_image[x][y];
            B = (Bx < 0 || By < 0)? 0 : pt_image[Bx][By];
            C = (Cx < 0 || Cy < 0)? 0 : pt_image[Cx][Cy];
            D = (Dx < 0 || Dy < 0)? 0 : pt_image[Dx][Dy];
            if (A == 0)
            {
                label_image[x][y] = 0;
            }else if ((B == 0) && (C == 0) && (D == 0))
            {
                label_image[x][y] = label;
                equiv_table[label] = label;
                label += 1;
            }else if( (D > 0))
            {
                label_image[x][y] = label_image[Dx][Dy];
                
            }else if((D == 0) && (B == 0) && (C > 0))
            {
                label_image[x][y] = label_image[Cx][Cy];
            }else if ((D == 0) && (C == 0) && (B > 0))
            {
                label_image[x][y] = label_image[Bx][By];
            }else if ((D == 0) && (B > 0) && (C > 0))
            {
                if (label_image[Cx][Cy] == label_image[Bx][By])
                {
                    label_image[x][y] = label_image[Bx][By];
                }else
                {
                    if (label_image[Cx][Cy] > label_image[Bx][By])
                    {
                        equiv_table[label_image[Cx][Cy]] = equiv_table[label_image[Bx][By]];
                    }else
                    {
                        equiv_table[label_image[Bx][By]] = equiv_table[label_image[Cx][Cy]];
                    }
                    label_image[x][y] = label_image[Bx][By];
                }                   
            }  
        }       
    }
    int num = label_image[0][0];
    for (int x = 0; x < height; x++)
    {
        for (int y = 0; y < width; y++)
        {
            label_image[x][y] = equiv_table[label_image[x][y]];
            num = label_image[x][y] > num ? label_image[x][y] : num;
        }        
    }
    int val, rept, max = 0;
    for (size_t i = 1; i < (num + 1); i++)
    {
        val = equiv_table[i];
        rept = 0;
        for (size_t j = 0; j < i; j++)
        {
            if (val == set[j])
            {
                rept = 1;
            }            
        }
        if (!rept)
        {
            max += 1;
            set[max] = val;
            
        }      
    }    
    *(num_of_objects) = max;
    return (label_image);
}
double** image_extract_getobj(short **label_image, double **pt_image, int height, int width, int label)
{
    double **nimage = image_create_image_channel(height,width);
    for (int x = 0; x < height; x++)
    {
        for (int y = 0; y < width; y++)
        {
            if (label_image[x][y] == label)
            {
                nimage[x][y] = pt_image[x][y];
            }
        }       
    }
    return (nimage);
}
void image_extract_save_objs(short **label_image, double **pt_image, int height, int width, short *set,int num_of_objects, char *generalname)
{
    double ***obj = NULL;
    char name[200] = {0};   
    for (size_t i = 1; i <= num_of_objects; i++)
    {
        obj = (double ***)malloc(sizeof(double **) * 1);
        obj[0] = image_extract_getobj(label_image,pt_image,height,width,set[i]);
        sprintf(name, "%d_%s",i,generalname);
        printf("saving object: %d\n", i);
        support_save_image(obj,width,height,1,name);
    }
    return;
}

void image_extract_objects(double ***image, int height, int width, int num_of_channels, char *argv)
{
    //segment image
    
    //smooth image
    double ***smoothed_image = image_spatial_smooth(image,height,width,num_of_channels,3,3,1);
    image_free_image(image,height,num_of_channels);
    double ***new_image = image_create_image(height,width,num_of_channels);
    for (size_t i = 0; i < num_of_channels; i++)
    {
        new_image[i] = image_seg_global_thres(smoothed_image[i],height,width,0,100);
    }
    int num_of_objects = 0;
    short *set = calloc(1000,sizeof(short));
    short **label_image = image_extract_label_regions(new_image[0],height,width,set,&num_of_objects);
    printf("Done Here\n");
    image_free_image(new_image,height,num_of_channels);
    image_extract_save_objs(label_image, smoothed_image[0],height,width,set, num_of_objects,argv);
    free(set);
    image_free_image(smoothed_image,height,num_of_channels);
    for (size_t i = 0; i < height; i++)
    {
       free(label_image[i]);
    }
    free(label_image);
    return;
}

double image_extract_getarea(double **binregion,int height, int width)
{
    double area = 0;
    for (int x = 0; x < height; x++)
    {
        for (int y = 0; y < width; y++)
        {
            if (binregion[x][y] > 0)
            {
                area += 1;
            }
            
            
        }       
    }
    return (area);
}
int timert = 0;
double image_extract_getperimeter(double **binregion, int height, int width)
{
    double B[3][3] = {{1,1,1},{1,1,1},{1,1,1}};
    double **boundary = image_boundary_morphology(binregion,height, width,3, 3, B);
    double perimeter = 0;
    //char name[200] = {0};  
   //sprintf(name, "boundary_%d_.png",timert);
    for (int x = 0; x < height; x++)
    {
        for (int y = 0; y < width; y++)
        {
            if (boundary[x][y] > 0)
            {
                perimeter += 1.0;
            }            
        }       
    }
    //double ***obj = (double ***)malloc(sizeof(double **) * 1);
    //obj[0] = boundary;
    //support_save_image(obj,width,height,1,name);
    //image_free_kernel(boundary,height);
    return (perimeter);
}
void image_extract_getcentroid(double **binregion, int height, int width, double *cx, double *cy, double area)
{
    double ybar = 0, xbar = 0;
    for (int x = 0; x < height; x++)
    {
        for (int y = 0; y < width; y++)
        {
            if (binregion[x][y] > 0)
            {
                ybar += y ;
                xbar += x ;
            }            
        }       
    }
    if (area > 0)
    {
        (*cx) = xbar / area;
        (*cy) = ybar / area;
    }else
    {
        (*cx) = 0;
        (*cy) = 0;
    }
        
}
void image_print_regoion_features(double **binregion, int height, int width)
{
    double area = image_extract_getarea(binregion,height,width);
    double perimeter = image_extract_getperimeter(binregion,height,width);
    double compactness = 999999999999999999, circularity = 9999999999999,effect_dia;
    if (area > 0)
    {
        compactness = (perimeter * perimeter) / area;
    }
    if (perimeter > 0)
    {
        circularity = (4 * 3.141592654 * area) / (perimeter * perimeter);
    }
    effect_dia = 2 * sqrt(area / 3.141592654);
    double convariance[2][2];
    image_extract_get_covariance_matrix(binregion,height,width,convariance);
    double eigen_one, eigentwo, eccentricity = 600, lambda;
    int eig = image_extract_get_eigen_values(binregion,height, width,convariance,&eigen_one,&eigentwo);
    if (eig)
    {
        lambda = eigen_one;
        if (eigen_one < eigentwo)
        {
            eigen_one = eigentwo;
            eigentwo = lambda;
        }
        if (eigen_one > 0)
        {
            lambda = eigentwo / eigen_one;
            lambda *= lambda;
            eccentricity = sqrt(1.0 - lambda);
        }
    }
    
    double cx, cy;
    image_extract_getcentroid(binregion,height,width,&cx,&cy,area);
    printf("centroid          : (%6.2f,%-6.2f)\n",cx,cy);
    printf("perimeter         : %10.5f\n",perimeter);
    printf("area              : %10.5f\n",area);
    printf("compactness       : %10.5f\n",compactness);
    printf("circularity       : %10.5f\n",circularity);
    printf("effective diameter: %10.5f\n",effect_dia);
    printf("convariance matrix: \n");
    for (size_t i = 0; i < 2; i++)
    {
        printf("                   ");
        for (size_t j = 0; j < 2; j++)
        {
            printf("%10.5f ",convariance[i][j]);
        }
    printf("\n");
    }
    printf("eigen values      :(%10.5f,%-10.5f)\n",eigen_one,eigentwo);
    printf("ecentricity       :%10.5f\n\n",eccentricity);
    return;
}
void image_extract_objinfo(double ***image, int height, int width, int num_of_channels)
{
    double ***smoothed_image = image_spatial_smooth(image,height,width,num_of_channels,3,3,1);
    image_free_image(image,height,num_of_channels);
    double ***new_image = image_create_image(height,width,num_of_channels);
    for (size_t i = 0; i < num_of_channels; i++)
    {
        new_image[i] = image_seg_global_thres(smoothed_image[i],height,width,0,100);
    }
    image_free_image(smoothed_image,height,num_of_channels);
    int num_of_objects = 0;    
    short *set = calloc(1000,sizeof(short));
    short **label_image = image_extract_label_regions(new_image[0],height,width,set,&num_of_objects);
    double **region = NULL;  
    for (size_t i = 1; i <= num_of_objects; i++)
    {
        region = image_extract_getobj(label_image,new_image[0],height,width,set[i]);
        printf("Region(%d) features:\n",i);
        image_print_regoion_features(region,height,width);  
        image_free_kernel(region,height);   
        timert +=1;
    }
    free(set);
    image_free_image(new_image,height,num_of_channels);
    for (size_t i = 0; i < height; i++)
    {
       free(label_image[i]);
    }
    free(label_image);
    return;
}
void image_extract_save_objinfo(double ***image, int height, int width, int num_of_channels,char *filename ,int label)
{
    double ***smoothed_image = image_spatial_smooth(image,height,width,num_of_channels,3,3,1);
    image_free_image(image,height,num_of_channels);
    double ***new_image = image_create_image(height,width,num_of_channels);
    for (size_t i = 0; i < num_of_channels; i++)
    {
        new_image[i] = image_seg_global_thres(smoothed_image[i],height,width,0,100);
    }
    image_free_image(smoothed_image,height,num_of_channels);
    char name[256] = {0};
    char *dot = strchr(filename,'.');
    int n = dot - filename;
    strncpy(name,filename,n);
    sprintf(name, "%s.txt",name);
    FILE *fp = fopen(name,"w");
    int num_of_objects = 0;    
    short *set = calloc(1000,sizeof(short));
    short **label_image = image_extract_label_regions(new_image[0],height,width,set,&num_of_objects);
    double **region = NULL;  
    for (size_t i = 1; i <= num_of_objects; i++)
    {
        region = image_extract_getobj(label_image,new_image[0],height,width,set[i]);
        image_extract_save_region_features(region,height,width,fp,label);
        image_free_kernel(region,height);
    }
    free(set);
    image_free_image(new_image,height,num_of_channels);
    for (size_t i = 0; i < height; i++)
    {
       free(label_image[i]);
    }
    free(label_image);
    return;
}
void image_extract_save_region_features(double **binregion, int height, int width,FILE *fp, int i)
{
    double area = image_extract_getarea(binregion,height,width);
    double perimeter = image_extract_getperimeter(binregion,height,width);
    double compactness = 0, circularity = 0;
    if (area > 0)
    {
        compactness = (perimeter * perimeter) / area;
    }
    if (perimeter > 0)
    {
        circularity = (4 * 3.141592654 * area) / (perimeter * perimeter);
    }
    fprintf(fp,"%f,%f,%d\n",compactness,circularity,i);
    return;
}

/**
 * @brief finds the eigen values of the covariance matrix of a region
*/
int image_extract_get_eigen_values(double **binregion, int height, int width, double (*CONV)[2],double *eigen_one, double *eigentwo)
{
    double a,b,c,d, under_root ,numer;
    b = CONV[0][0] + CONV[1][1];
    c = 4 * ((CONV[0][0] * CONV[1][1]) - (CONV[0][1] * CONV[0][1]));
    under_root = b * b - c;
    if (under_root < 0)
    {
        return (0);
    } 
    numer = sqrt(under_root);
    (*eigen_one) = (b + numer) / 2;
    (*eigentwo) = (b - numer) / 2;
    return (1);
    
}

/**
 * @brief finds the convariance matrix of a region
*/
void image_extract_get_covariance_matrix(double **binregion, int height, int width, double (*CONV)[2])
{
    double mean[2] = {0};
    double convariance[2][2] = {0};
    CONV[0][0] = 0; CONV[0][1] = 0; CONV[1][0] = 0; CONV[1][1] = 0;
    int num = 0;
    //compute the mean 
    for (int x = 0; x < height; x++)
    {
        for (int y = 0; y < width; y++)
        {
            if (binregion[x][y] > 0)
            {
                mean[0] += x;
                mean[1] += y;
                num++;
            }            
        }       
    }
    if (num > 0)
    {
        mean[0] /= num;
        mean[1] /= num;
    }
    double xdiff,ydiff;
    for (int x = 0; x < height; x++)
    {
        xdiff = x - mean[0];
        for (int y = 0; y < width; y++)
        {
            if (binregion[x][y] > 0)
            {
                ydiff = y - mean[1];
                CONV[0][0] += xdiff * xdiff;
                CONV[0][1] += xdiff * ydiff;
                CONV[1][0] += xdiff * ydiff;
                CONV[1][1] += ydiff * ydiff;
            }            
        }       
    }
    double K = num - 1;
    if (K > 1)
    {
        CONV[0][0] /= K;
        CONV[0][1] /= K;
        CONV[1][0] /= K;
        CONV[1][1] /= K;
    }
    return;
}

/**
 * @brief harris stephen corner detector
*/
double **image_extract_harris_corner_detector(double **img, int height, int width,int m, int n, double k, double thres)
{  
    double **f_x = image_create_image_channel(height,width);
    double **f_y = image_create_image_channel(height,width);
    double **f_c = image_create_image_channel(height,width);
    //compute the derivatives
    int xup, xdown, yleft, yright, tx, ty;
    double xdev, ydev;
    for (size_t x = 1; x < height - 1; x++)
    {
        xup = x - 1;        
        xdown = x + 1;
        for (size_t y = 1; y < width - 1; y++)
        {
            yleft = y - 1;
            yright = y + 1;
            xdev = img[xdown][y] - img[xup][y];
            ydev = img[x][yright] - img[x][yleft];
            f_x[x][y] = xdev;
            f_y[x][y] = ydev;
        }      
    }
    int i, j;
    double a,b,det, trace, r;
    double M[2][2] = {0};
    a = m / 2;
    b = n / 2;
    for (size_t x = 1; x < height - 1; x++)
    {
        for (size_t y = 1; y < width - 1; y++)
        {
            M[0][0] = M[0][1] = M[1][0] = M[1][1] = 0.0;
            for (int s = -a; s <= a; s++)
            {
                i = x - s;              
                for (int t = -b; t <= b; t++)
                {
                    j = y - t;
                    xdev = f_x[i][j];
                    ydev = f_y[i][j];
                    M[0][0] += xdev * xdev;
                    M[0][1] += xdev * ydev;
                    M[1][0] += xdev * ydev;
                    M[1][1] += ydev * ydev;
                } 
            }            
            det = (M[0][0] * M[1][1]) - (M[0][1] * M[0][1]);
            trace = M[0][0] + M[1][1];            
            r = det - k * trace * trace;
            if (r > thres)
            {
                f_c[x][y] = 255.0;
            }            
        }      
    }
    image_free_kernel(f_x,height);
    image_free_kernel(f_y,height);
    return (f_c);
}

/**
 * @brief implementation of SIFT algorithm
*/
void image_extract_gaussians(double **image, int height , int width, double std, int s , char (*names)[15])
{
    int kernel_dim;
    double *gaus_kernel = NULL;
    int num_smoothed_imgs = s + 3;
    double k = pow(2, 1.0/ ((double)s));
    printf("k = %f\n",k);
    double p;
    double **smoothed = NULL;
    double ***nimage = NULL;
    FILE *fp;
    for (size_t i = 0; i < num_smoothed_imgs; i++)
    {
        sprintf(&(names[i][0]), "smoothed_%d.png",i);
        p = pow(k, i);
        printf("std = %10.5f\n",p * std);
        //gaus_kernel = image_extract_unsep_gauss(&kernel_dim,p * std);
        //smoothed = image_spatial_filter(image,gaus_kernel,height,width,1,kernel_dim,kernel_dim);
        gaus_kernel = image_extract_sep_gauss(&kernel_dim, p * std);
        smoothed = image_extract_gaus_smooth(image,gaus_kernel,gaus_kernel,height, width,kernel_dim, kernel_dim);
        free(gaus_kernel);
        //image_free_kernel(gaus_kernel,kernel_dim);
        fp = fopen(&(names[i][0]), "wb");
        if (fp == NULL)
        {
            printf("Failure\n");
        }
        //support_save_image(smoothed[0],width,height,1,&(names[i][0]));
        support_temporal_save_image(smoothed,height,width,fp);
        image_free_kernel(smoothed,height);
        //image_free_image(smoothed,height,1);
    }
    return;   
}
void image_extract_gausdifference(int *height , int *width, char (*names)[15], char (*gausdiff)[15], int num_images)
{
    //read image 1;
    double **image_one, **image_two, **dog, diff;
    double ***nimage = NULL;
    FILE *fp = NULL;
    fp = fopen(&(names[0][0]),"rb");
    image_one = support_get_temporal_saved_image(height, width,fp);
    int h = (*height);
    int w = (*width);
    printf("h = %3d : w = %3d\n",h,w);
    for (size_t i = 1; i < num_images; i++)
    {
        fp = fopen(&(names[i][0]),"rb");
        image_two = support_get_temporal_saved_image(height,width,fp);
        nimage = malloc(sizeof(double **));
        nimage[0] = image_one;
        sprintf(&(gausdiff[i][0]),"dog_%d.png",i - 1);
        dog = image_create_image_channel(h,w);
        // compute the difference
        for (size_t x = 0; x < h; x++)
        {
            for (size_t y = 0; y < w; y++)
            {
                diff = image_two[x][y] - image_one[x][y];
                diff = diff < 0 ? -diff: diff;
                dog[x][y] = diff;           
            }
        }
    
           /*for (size_t y = 0; y < w; y++)
            {
                printf("%3f ",dog[h/2][y]);          
            }
            printf("\n");*/ 
        //save the dog to file.
        //fp = fopen(gausdiff[i],"wb");
        nimage = malloc(sizeof(double **));
        nimage[0] = dog;
        image_correct_dynamic_range(nimage,h,w,1);
        support_save_image(nimage,w,h,1,&(gausdiff[i][0]));
        //support_temporal_save_image(dog,height,width,fp);
        //image_free_kernel(dog,height);
        image_free_kernel(image_one,h);
        image_one = image_two;
    }
    return;   
}
/**
 * @brief performs gaussian smoothing of an image
*/
double **image_extract_gaus_smooth(double **image,double *kernel_row,double *kernel_col,size_t image_height, size_t image_width,int m , int n)
{
    double **new_image = image_create_image_channel(image_height,image_width);
    int a = (m - 1) / 2;
    int b =(n - 1) / 2;
    double gxy;
    int i;
    int j;
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
            new_image[x][y]  = gxy; 
        }      
    }
    double **new_image_2 = image_create_image_channel(image_height,image_width);
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
          gxy += kernel_col[t + b] * new_image[j][y];       
        }
        new_image_2[x][y]  = gxy; 
      }      
    } 
  image_free_kernel(new_image,image_height);
  return (new_image_2);
}

double *image_extract_sep_gauss(int *m,double std)
{
  double divisor = 2 * 3.141592654 * std *std;
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
    gaus[t + c] = exp( ((double)( t * t) )/ (mul)) / divisor;
    sum += gaus[t + c];
  }
  //printf("m = %d\n",*m);
  for (int i = 0; i < max; i++)
  {
    gaus[i] /= sum;
  }
  //printf("sum = %f\n",sum);
  double **kernel = image_reconstruct_kernel(gaus,gaus,max,max);
  printf("kernel: \n");
  sum = 0;
  for (size_t i = 0; i < max; i++)
  {
    for (size_t j = 0; j < max; j++)
    {
        printf("%10.5f ",kernel[i][j]);
        sum += kernel[i][j];
    }
    printf("\n");
  }
  printf("sum = %f\n",sum);
  image_free_kernel(kernel,max);
  return (gaus);

}