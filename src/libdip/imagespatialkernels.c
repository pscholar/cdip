double *image_create_box_kernel_sep(int m)
{
  double *box = (double *)malloc(sizeof(double) * m);
  double divisor = m;
  double coefficient = 1.0F / divisor;
  for (size_t i = 0; i < m; i++)
  {
    box[i]= coefficient;
  }
  return (box);
}
double **image_create_box_kernel_un(int m, int n)
{
  double **box = (double **)malloc(sizeof(double *) * m);
  double divisor = m * n;
  double coefficient = 1.0F / divisor;
  for (size_t i = 0; i < m; i++)
  {
    box[i]= (double *)malloc(sizeof(double) * n);
    for (size_t j = 0; j < n; j++)
    {
      box[i][j] = coefficient;
    }
    
  }
  return (box);
}
double **image_create_gaussian_kernel_un(int m,int n, double k)
{
  double sd = ( ((double)m )/ 6.0F) * ( ((double)m )/ 6.0F) * 2;
  double **gauss = (double **)malloc(sizeof(double *) * m);
  double sum = 0.0F;
  double value;
  double x;
  int a = (m - 1) / 2;
  int b = (n - 1) / 2;
  for (int s = -a; s <= a; s++)
  {
    gauss[s + a]= (double *)malloc(sizeof(double) * n);
    x = ((double) (s * s));
    for (int t = -b;  t<= b; t++)
    {
      value = k * exp(-(x + ((double) (t * t))) / sd);
      gauss[s + a][ t + b] = value;
      sum += value;
    }
    
  }
  for (size_t i = 0; i < m; i++)
  {
    for (size_t j = 0; j < n; j++)
    {
      gauss[i][j] = gauss[i][j] / sum;
    }

  }
  return (gauss);

}
double *image_create_gaussian_kernel_sep(int m,double k)
{
  double sd = ( ((double)m )/ 6.0F) * ( ((double)m )/ 6.0F) * 2;
  sd = 2.0f;
  double *gauss = (double *)malloc(sizeof(double) * m);
  double sum = 0.0F;
  double value;
  int a = (m - 1) / 2;
  for (int s = -a; s <= a; s++)
  {
      value = k * exp(-(((double) (s * s))) / sd);
      gauss[s + a] = value;
      sum += value;
    
  }
  for (size_t i = 0; i < m; i++)
  {
    gauss[i] = gauss[i] / sum;
      
  }
  return (gauss);

}

void image_free_kernel(double **kernel,int m)
{
  for (size_t i = 0; i < m; i++)
  {
    free(kernel[i]);
  }
  free(kernel);
  return; 
}

double ** image_create_sharpfilter(double *w, double *v,int m, int n)
{
  // create an impulse of size m x n.
  double **impulse = image_create_impulse(m, n);
  printf("Impulse: \n");
  mathlib_matrixlib_printarray(impulse,m,n);
  double **lp = mathlib_matrixlib_create_zero(m,n);
  mathlib_matrixlib_vector_multiplication(w,v,lp,m,n);
  printf("low pass filter: \n");
  mathlib_matrixlib_printarray(lp,m,n);
  //create a high pass filter.
  double **hp = mathlib_matrixlib_create_zero(m,n);
  mathlib_matrixlib_subtraction(impulse,lp,hp,m,n);
  printf("high pass filter: \n");
  mathlib_matrixlib_printarray(hp,m,n);
  image_free_kernel(lp,m);
  return (hp);

}
double **image_create_bandreject_filter(double **lp1, double** lp2, int m, int n)
{
  double **impulse = image_create_impulse(m, n);
  printf("Impulse: \n");
  mathlib_matrixlib_printarray(impulse,m,n);
  double **br = mathlib_matrixlib_create_zero(m,n);
  mathlib_matrixlib_subtraction(impulse,lp2,br,m,n);
  mathlib_matrixlib_addition(lp1,br,br,m,n);
  printf("band reject filter: \n");
  mathlib_matrixlib_printarray(br,m,n);
  return (br);

}
double **image_create_bandpass_filter(double **br, int m, int n)
{
  double **impulse = image_create_impulse(m, n);
  printf("Impulse: \n");
  mathlib_matrixlib_printarray(impulse,m,n);
  double **bp = mathlib_matrixlib_create_zero(m,n);
  mathlib_matrixlib_subtraction(impulse,br,bp,m,n);
  printf("band pass filter: \n");
  mathlib_matrixlib_printarray(bp,m,n);
  return (bp);
}
double **image_reconstruct_kernel(double *w, double *v,int m, int n)
{
  double **lp = mathlib_matrixlib_create_zero(m,n);
  mathlib_matrixlib_vector_multiplication(w,v,lp,m,n);
  return (lp);
}
void image_split_kernel(double **kernel,double *w, double *v , int m, int n)
{
  //search for non-zero coefficient in the array.
  double coefficient = 1.0F;
  double zero = 0.0F;
  int x , y;
  int isFound = 0;
  for (x = 0; x < m; x++)
  {
    for (y = 0; y < n; y++)
    {
      if ( (kernel[x][y] > zero) || (kernel[x][y] < zero))
      {
        coefficient = kernel[x][y];
        isFound = 1;
        break;
      }
    }
    if (isFound)
    {
      break;
    }
    
  }
  //get the row.
  if (isFound)
  {
    printf("column vector: \n");
    for (size_t j = 0; j < m; j++)
    {
      w[j] = kernel[j][y];
      printf("%8.5f\n",w[j]);
    }
    printf("row vector:\n");
    for (size_t i = 0; i < n; i++)
    {
      v[i] = kernel[x][i] / coefficient;
      printf("%8.5f ",v[i]);
    }
    printf("\n");
  }
  double **k = image_reconstruct_kernel(w,v,m,n);
  printf("reconstructed kernel: \n");
  mathlib_matrixlib_printarray(k,m,n);
  image_free_kernel(k,m);
  return;

}
double **image_create_impulse(int m, int n)
{
  double **impulse = (double **)malloc(sizeof(double *)*m);
  for (size_t i = 0; i < m; i++)
  {
    impulse[i] = (double *)calloc(n, sizeof(double));
  }
  int cx = m / 2;
  int cy = n / 2;
  impulse[cx][cy] = 1;
  return impulse;
}
