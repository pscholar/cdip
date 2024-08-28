#include "matrixlib.h"


/**
 * @brief prints the elements of a 2D array.
 * @param M number of rows.
 * @param N number of rows.
*/
void mathlib_matrixlib_printarray(double **A, int M, int N)
{
  for (size_t i = 0; i < M; i++)
  {
    for (size_t j = 0; j < N; j++)
    {
      printf("%8.5f ", A[i][j]);
    }
    printf("\n");
  }
  return;
}


/**
 * Copies elements of one matrix into another of the same dimensions.
 * @param src: the source matrix.
 * @param dest: the destination matrix.
 * @param M: number of rows.
 * @param N: number of columns.
*/
void mathlib_matrixlib_copy(double **src, double **dest, size_t M, size_t N)
{
  for (size_t i = 0; i < M; i++)
  {
    for (size_t j = 0; j < N; j++)
    {
      dest[i][j] = src[i][j];
    }
  }
  return;
}

/**
 * Copies pointers to rows of a matrix into an array of pointers.
 * @param src: the source matrix.
 * @param dest: the destination matrix.
 * @param M: number of rows.
 * @param N: number of columns.
*/
void mathlib_matrixlib_copy_pointers(double *src, double **dest, size_t M, size_t N)
{
  size_t row;
  printf("printing the pointers\n");
  for (size_t i = 0; i < M; i++)
  {
    row = i * N;
    dest[i] = (src + row);
    
  }  
  return;
}

/**
 * Does: transposes a matrix
 * @param A: pointer to 2D array
 * @param C: pointer to 2D array
 * @param M: number of rows.
 * @param N: number of columns
 * @implements: C = A'
*/
void mathlib_matrixlib_transpose(double *A, double *B, size_t M, size_t N)
{
  size_t row_A;
  size_t row_B;
  for (size_t i = 0; i < M; i++)
  {
    row_A = i * N;  
    for (size_t j = 0; j < N; j++)
    {
      row_B = j * M;
      *(B + row_B + i) = (*(A + row_A + j));
    }
    
  }
  return;
}

/**
 * Does: element-wise matrix addition
 * @param A: pointer to 2D array
 * @param B: pointer to 2D array
 * @param C: pointer to 2D array
 * @param M: number of rows.
 * @param N: number of columns
 * @implements: C = A + B
*/
void mathlib_matrixlib_addition( double **A, double **B, double **C, size_t M, size_t N )
{
  for (size_t i = 0; i < M; i++)
  {

    for (size_t j = 0; j < N; j++)
    {
      C[i][j] = A[i][j] + B[i][j];
    }
    
  }
  return;

}

/**
 * @brief: Does element-wise matrix addition, with scalar multiplcation provision;
 * @param A: pointer to 2D array
 * @param B: pointer to 2D array
 * @param C: pointer to 2D array
 * @param k_A: scalar to multiply with A.
 * @param k_B: scalar to multiply with B.
 * @param M: number of rows.
 * @param N: number of columns
 * @implements: C = k_A * A + k_B * B
*/
void mathlib_matrixlib_addition_scalar( double **A, double **B, double **C, double k_A, double k_B,size_t M, size_t N )
{
  for (size_t i = 0; i < M; i++)
  {
    for (size_t j = 0; j < N; j++)
    {
      C[i][j] =  k_A * A[i][j] + k_B * B[i][j];
    }
    
  }
  return;

}

/**
 * Does: element-wise matrix subtraction
 * @param A: pointer to 2D array
 * @param B: pointer to 2D array
 * @param C: pointer to 2D array
 * @param M: number of rows.
 * @param N: number of columns
 * @implements: C = A - B
*/
void mathlib_matrixlib_subtraction( double **A, double **B, double **C, size_t M, size_t N )
{
  mathlib_matrixlib_addition_scalar(A,B,C,1.0F, -1.0F, M,N);
  return;

}

/**
 * @brief: Does element-wise matrix subtraction, with scalar multiplcation provision;
 * @param A: pointer to 2D array
 * @param B: pointer to 2D array
 * @param C: pointer to 2D array
 * @param k_A: scalar to multiply with A.
 * @param k_B: scalar to multiply with B.
 * @param M: number of rows.
 * @param N: number of columns
 * @implements: C = k_A * A - k_B * B
*/
void mathlib_matrixlib_subtraction_scalar( double **A, double **B, double **C, double k_A, double k_B,size_t M, size_t N )
{
  mathlib_matrixlib_addition_scalar(A,B,C,k_A, -k_B, M,N);
  return;

}

/**
 * @brief: Does element-wise scalar multiplication;
 * @param A: pointer to 2D array
 * @param B: pointer to 2D array
 * @param k: scalar.
 * @param M: number of rows.
 * @param N: number of columns
 * @implements: B =  kA
*/
void mathlib_matrixlib_scalar_multplication( double **A, double **B, double k,size_t M, size_t N)
{
  for (size_t i = 0; i < M; i++)
  {
    for (size_t j = 0; j < N; j++)
    {
      B[i][j] = k * A[i][j];
    }
    
  }

  return;
}
/**
 * @brief performs multiplication of two comformable vectors.
 * @param w column vector of dimensions m x 1.
 * @param v row vector of dimensions n x 1.
 * @param m number of rows in w.
 * @param n number of columns in v.
*/
void mathlib_matrixlib_vector_multiplication(double *w, double *v,double **r ,int m, int n)
{
   for (size_t k = 0; k < n; k++)
  {
    for (size_t i = 0; i < m; i++)
    {
      r[i][k] = w[i] * v[k];
    }
  
  }
  return;
}

/**
 * @brief: Does element-wise matrix multiplication;
 * @param A: pointer to 2D array
 * @param B: pointer to 2D array
 * @param C: pointer to 2D array
 * @param M: number of rows.
 * @param N: number of columns
 * @implements: C =  A .* B
*/
void mathlib_matrixlib_elem_multplication( double **A, double **B, double **C,size_t M, size_t N)
{
  for (size_t i = 0; i < M; i++)
  {
    for (size_t j = 0; j < N; j++)
    {
      C[i][j] = A[i][j]* B[i][j];
    }
    
  }
  return;
}

/**
 * @brief: Does element-wise matrix multiplication, with scalar multiplcation provision;
 * @param A: pointer to 2D array
 * @param B: pointer to 2D array
 * @param C: pointer to 2D array
 * @param k_A: scalar to multiply with A.
 * @param k_B: scalar to multiply with B.
 * @param M: number of rows.
 * @param N: number of columns
 * @implements: C = k_A * A * k_B * B
*/
void mathlib_matrixlib_elem_multplication_scalar( double **A, double **B, double **C, double k_A, double k_B,size_t M, size_t N)
{
  double scalar_pdct = k_A * k_B;
  for (size_t i = 0; i < M; i++)
  {
    for (size_t j = 0; j < N; j++)
    {
      C[i][j] = scalar_pdct*A[i][j]* B[i][j];
    }
  } 
  return;
}

/**
 * Does: does normal matrix multiplication
 * @param A: 2D matrix of dimensions M x N
 * @param B: 2D matrix of dimensions N x R
 * @param C: destination 2D matrix of dimensions M * R
 * @param M: number of rows in matrix A.
 * @param N: number of columns in A same as number of rows in B.
 * @param R: number of columns in B.
*/
void mathlib_matrixlib_multiplication(double **A, double **B, double **C,size_t M, size_t N, size_t R)
{

  for (size_t k = 0; k < R; k++)
  {
    for (size_t i = 0; i < M; i++)
    {
      C[i][k] = 0;
      for (size_t j = 0; j < N; j++)
      {
        C[i][k] += A[i][j] * B[j][k];
      }
    }
  
  }
  return;

}

/**
 * Does: does normal matrix multiplication, with scalar provision
 * @param A: 2D matrix of dimensions M x N
 * @param B: 2D matrix of dimensions N x R
 * @param C: destination 2D matrix of dimensions M * R
 * @param k_A: scalar for matrix A
 * @param k_B: scalar for matrix B.
 * @param M: number of rows in matrix A.
 * @param N: number of columns in A same as number of rows in B.
 * @param R: number of columns in B.
*/
void mathlib_matrixlib_multiplication_scalar(double **A, double **B, double **C,double k_A, double k_B, size_t M, size_t N, size_t R)
{
  double scalar = k_A * k_B;
  for (size_t k = 0; k < R; k++)
  {
    for (size_t i = 0; i < M; i++)
    {
      C[i][k] = 0;
      for (size_t j = 0; j < N; j++)
      {
        C[i][k] += scalar *  A[i][j] * B[j][k];
      }
    }
  
  }
  return;

}

/**
 * @brief: Finds the trace of a matrix;
 * @param A: matrix of dimensions M x M
 * @param M: dimension of the matrix.
 * @implements: trace(A)
*/
double mathlib_matrixlib_trace(double **A, size_t M)
{
  double trace = 0;
  // sum along the diagonals.
  for (size_t i = 0; i < M; i++)
  {
    trace += A[i][i];
  }
  
  return (trace);
}

/**
 * @brief: generates an identity matrix.
 * @param M: dimension of the identity matrix.
 * @returns: pointer to the matrix.
*/
double** mathlib_matrixlib_create_unit(size_t M)
{
  double **identity = (double **)malloc(sizeof(double *) * M);
  for (size_t i = 0; i < M; i++)
  {
    identity[i] = (double *)calloc(M, sizeof(double));
    identity[i][i] = 1.0F;
  }
  return (identity);
}
/**
 * @brief generates a zero matrix.
 * @param m number of rows.
 * @param n number of columns
 * @returns: pointer to the matrix.
*/
double** mathlib_matrixlib_create_zero(int m, int n)
{
  double **zero = (double **)malloc(sizeof(double *) * m);
  for (size_t i = 0; i < m; i++)
  {
    zero[i] = (double *)calloc(n, sizeof(double));
  }
  return (zero);

}
/**
 * Swaps a row whose leading element in a particular column is 0 with a row with a non-zero leading element if possible.
 * @param  A: array of pointers to rows.
 * @param M: dimension of A.
 * @param index: row and column number we are interested in.
 * @param mul_factor: factor to be multiplied by a negative in case of a swap.
*/
void mathlib_matrixlib_swaprow(double **A,size_t M,size_t index, double *mul_factor)
{
  double *temp;  
  if (A[index][index] == 0.0F )
  {
    for (int k = index; k < M; k++)
    {
      if (A[k][index] != 0.0F)
      {
        /*Make the swap*/
        temp = A[index];
        A[index] = A[k];
        A[k] = temp;
        (*(mul_factor)) *= -1.0F;
        //printf("swapped\n");
        break;
      }        
    }
  }
  return;
    
}

/**
 * Reduces a square matrix to an echelon form / triangular matrix.
 * @param A: pointer to the rows of the matrix.
 * @param M: number of rows in the matrix.
 * @param N: number of columns in the matrix.
 * @returns: the multiplicative number.
*/
double  mathlib_matrixlib_row_reduction(double **A, size_t M,size_t N)
{
  double mul_factor = 1.0F;
  for (size_t i = 0; i < M - 1; i++)
  {
    /*for each leading element of a row*/
    mathlib_matrixlib_swaprow(A,M,i,&mul_factor);
    if (A[i][i] != 0.0F)
    {
      /*reduce all element below this to zero*/
      //printf("1.A[%ld][%ld] = %8.2f\n",i,i,A[i][i]);
      for (size_t j = (i + 1); j < M; j++)
      {
        //printf("2.A[%ld][%ld] = %8.2f\n",j,i,A[j][i]);
        if (A[j][i] != 0.0F)
        {
          double multiplier = ((A[j][i])/(A[i][i])*(-1.0F));
          //printf("Multiplier: %8.2f\n",multiplier);
          for (size_t k = i; k < N; k++)
          {
            A[j][k] = (A[i][k]) * multiplier + A[j][k];
            //printf("3.A[%ld][%ld] = %8.2f\n",j,k,A[j][k]);
          } 
        }
      }
    }   
  }

  return (mul_factor);
}

/**
 * @brief: computes the determinant of a matrix
 * @param A: square matrix of order M.
 * @param M: the order of the matrix.
 * @returns: the determinat of the matrix.
*/

double mathlib_matrixlib_det(double **A, size_t M)
{
  double mult = mathlib_matrixlib_row_reduction(A,M,M);
  double det = 1.0F;
  for (size_t i = 0; i < M; i++)
  {
    for (size_t j = 0; j < M; j++)
    {
      if(A[i][j] == 0.0F)
      {
        return (0.0F);
      }
      det *= A[i][j];
    }
  }
  
  return (det * mult);
}

/**
 * swaps rows of two augumented matrices laid side by side.
 * @param A: array of pointers to rows of a matrix of order M.
 * @param B: array of pointers to rows of a matrix of order M.
 * @param M: order of the matrices.
 * @param index: pivot under interest. 
*/
void mathlib_matrixlib_native_swaprow(double **A,double **B,size_t M,size_t index)
{
    double *temp;  
  if (A[index][index] == 0.0F )
  {
    for (int k = index; k < M; k++)
    {
      if (A[k][index] != 0.0F)
      {
        /*Make the swap*/
        temp = A[index];
        A[index] = A[k];
        A[k] = temp;
        temp = B[index];
        B[index] = B[k];
        B[k] = temp;
        //printf("swapped\n");
        break;
      }        
    }
  }
  return;
}

/**
 * reduces two augumented matrices to gauss-jordan form.
 * @param A: array of pointers to rows of a matrix of order M.
 * @param B: array of pointers to rows of a matrix of order M.
 * @param M: order of the matrices.
 * @param N: N = M. 
*/
void mathlib_matrixlib_gaussreduction(double **A,double **B,size_t M,size_t N)
{
  for (size_t i = 0; i < M; i++)
  {
    /*for each leading element of a row*/
    mathlib_matrixlib_native_swaprow(A,B,M,i);
    if (A[i][i] != 0.0F)
    {
      /*reduce all element below this to zero*/
      for (size_t j = 0; j < M; j++)
      {
        if (j != i)
        {
          if (A[j][i] != 0.0F)
          {
            double multiplier = ((A[j][i])/(A[i][i])*(-1.0F));
            for (size_t k = 0; k < N; k++)
            {
              A[j][k] = ((A[i][k]) * multiplier + A[j][k]);
              B[j][k] = (B[i][k]) * multiplier + B[j][k];
            } 
          }
        }

      }
    }   
  }
  for (size_t i = 0; i < M; i++)
  {
    if (A[i][i] == 0.0F)
    {
      printf("Non-Invertible");
      return;
    }
    for (size_t j = 0; j < M; j++)
    {
      B[i][j] /= A[i][i];
    }
    A[i][i] = 1.0F;   
  }
  
}

/**
 * calculates the inverse of a square matrix.
 * @param A: matrix of order M.
 * @param M: order of the matrix.
*/
double** mathlib_matrixlib_inverse(double **A, size_t M)
{

  //create an identity matrix.
  double **identity = mathlib_matrixlib_create_unit(M);
  mathlib_matrixlib_gaussreduction(A,identity,M,M);
  mathlib_matrixlib_copy(identity,A,M,M);
  return (identity);
}

/**
 * computes the rank of a matrix.
 * @param A: matrix of order M.
 * @param M: the order of the matrix.
*/
size_t mathlib_matrix_rank(double **A, size_t M)
{
  double **dup = mathlib_matrixlib_create_zero(M,M);
  //copy elements of A into its duplicate.
  mathlib_matrixlib_copy(A, dup, M, M);
  //reduce matrix to a diagonal matrix.
  mathlib_matrixlib_row_reduction(dup,M,M);
  //count all non-zero rows.
  size_t j = 0;
  size_t rank = 0;
  for (size_t i = 0; i < M; i++)
  {
    for (; j < M; j++)
    {
      // find the first non-zero element if this row.
      if (dup[i][j] != 0.0F)
      {
        rank += 1;
        j++;
        break;
      }
    }    
  }
  return (rank);
}
/**
 * The succedding versions of methods that allow specification of upper and lower limits when dealing with matrices
*/
/**
 * @brief: performs element-wise addition or subtraction.
 * @param A: matrix 
*/
void mathlib_matrix_linear_comb(double **A, double **B, double **C, double k_A, double k_B,size_t start_i,size_t start_j,size_t end_i,size_t end_j)
{
  for (size_t i = start_i; i < end_i; i++)
  {
    for (size_t j = start_j; j < end_j; j++)
    {
      C[i][j] =  k_A * A[i][j] + k_B * B[i][j];
    }
    
  }
}


/**
 * @brief stores a n-dimensional matrix in a file
 * @param matrix pointer to the matrix
 * @param rows number of rows in the matrix
 * @param cols number of columns in the matrix
 * @param n dimension of matrix
 * @param filename pointer to the file name
 * @returns pointer to a file
*/
FILE* mathlib_store_matrix(double ***matrix, size_t rows, size_t cols, int n, char *filename)
{
  FILE *fp;
  uint16_t header[] = {(uint16_t)rows,(uint16_t)cols,(uint16_t)n};
  size_t rowsize = cols * sizeof(double);
  if((fp = fopen(filename,"wb")) != NULL)
  {
    fwrite(header,(sizeof(header) / sizeof(header[0])) * sizeof(uint16_t),1,fp);
    for (int i = 0; i < n; i++)
    {
      for (size_t x = 0; x < rows; x++)
      {
        fwrite(matrix[i][x],rowsize,1,fp);
      }      
    }    
  }
  return fp;
}

/**
 * @brief extracts a channel of an n-matrix from a file
 * @param  rows address to rows variable
 * @param cols address to cols  varaible
 * @param n address to dimensionality variable
 * @param fp pointer to a file
 * @return a pointer to the nth channel.
*/
double** mathlib_retrieve_matrix_channel(size_t * rows, size_t *cols, int *n ,FILE *fp)
{
  uint16_t header[3];
  static int ip = 0;
  int rowsize;
  double **ptr = NULL;
  if (fp == NULL)
  {
    return ptr;
  }
  if (ip == 0)
  {
    fread(header,(sizeof(header) / sizeof(header[0])) * sizeof(uint16_t),1,fp);
    *rows = (size_t)header[0];
    *cols = (size_t)header[1];
    *n = (size_t)header[2];
    rowsize  = (*cols) * sizeof(double);
  }  
  
  if(ip < (*n))
  {
    ptr = (double **)malloc(sizeof(double *) * (*rows));
    for (size_t i = 0; i < (*rows); i++)
    {
      ptr[i] = (double *)malloc(sizeof(double) * (*cols));
    }
    for (size_t x = 0; x < (*rows); x++)
    {
      fread(ptr[x],rowsize,1,fp);
    }
    ip += 1;    
  }  
  return (ptr);
}

/**
 * @brief reads a matrix from a file
 * @param  rows address to rows variable
 * @param cols address to cols  varaible
 * @param n address to dimensionality variable
 * @param fp pointer to a file
 * @return a pointer to the matrix.
*/
double*** mathlib_matrix_retrieve_matrix(size_t * rows, size_t *cols, int *n ,FILE *fp)
{
  double ***image = NULL;
  double **ptr;
  int i = 0;
  if ((fp != NULL) && ((ptr = mathlib_retrieve_matrix_channel(rows,cols,n,fp) ) != NULL) )
  {
    image = (double ***)malloc(sizeof(double **)* (*n));
    image[i] = ptr;
    i++;
    while((ptr = mathlib_retrieve_matrix_channel(rows,cols,n,fp) ) != NULL)
    {
      image[i] = ptr;
      i++;
    }
  }
  return (image);
}
/*int main(void)
{
  double A[2][3] = {{1,2,3},{-1,0,2}};
  double B[2][3] = {{-1,5,-2},{2,2,-1}};
  printf("A: \n");
  printArray(&A[0][0],2,3);
  printf("B: \n");
  printArray(&B[0][0],2,3);
  double C[2][3];
  mathlib_matrixlib_addition(&A[0][0],&B[0][0],&C[0][0],2,3);
  printf("A + B; \n");
  printArray(&C[0][0],2,3);
  mathlib_matrixlib_addition_scalar(&A[0][0],&B[0][0],&C[0][0],1.0F,2.0F,2,3);
  printf("k1A + k2B; \n");
  printArray(&C[0][0],2,3);
  double AT[3][2];
  mathlib_matrixlib_transpose(&A[0][0],&AT[0][0], 2, 3);  
  printf("A' ;\n");
  printArray(&AT[0][0],3,2);
  double M[2][2] = {{4,-1},{0,2}};
  double N[2][3] = {{1,4,2},{3,1,5}};
  double E[2][3];
  mathlib_matrixlib_multiplication(&M[0][0], &N[0][0], &E[0][0],2,2,3);
  printf ("M x N: \n");
  printArray(&E[0][0],2,3);
  double T[2][2] = {{5,1},{1,26}};
  double trace = mathlib_matrixlib_trace(&T[0][0],2);
  printf("trace(T): %5.2f\n",trace);
  double *identity = mathlib_matrixlib_create_unit(3);
  printf("3 x 3 identity matrix: \n");
  printArray(identity,3,3);
  double F[3][3] = {{0,0,0},{2,3,0},{-1,-2,4}};
  printf("F:\n");
  printArray(&F[0][0],3,3);
  double det = mathlib_matrixlib_det(&F[0][0],3);
  printf("Det = %8.2f\n", det);
  double G[3][3] = {{3,2,1},{0,1,1},{-1,1,0}};
  mathlib_matrixlib_inverse(&G[0][0],3);
  return (0);
}
*/