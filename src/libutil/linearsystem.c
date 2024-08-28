#include <stdio.h>
#include <stdlib.h>
void swapRows(double *array[],int,int);
int rowReduction(double *array[],int,int );
void gaussJordan(double *echelon[],int,int);
double* particularSoln(double *reducedEchelon[], int row,int column);
void swapRows(double *array[],int row,int column)
{
     /**
      * Does: swaps row so that for each row, there is a leading variable that is non-zero
      * Recieves:
      * @array: an array of pointers, where each pointer points to a row of an augmented matrix.
      * @row: the number of rows in the augumented matrix.
      * @column: the number of columns in the coefficient matrix.
     */
     double *temp;
    for (int m = 0; m < (row - 1); m++)
     {
          
          if (array[m][m] == 0 )
          {
               /*scan through column n, and make a swap if there is a non-zero row*/
               for (int i = m; i < row; i++)
               {
                    if (array[i][m] > 0)
                    {
                         /*Make the swap*/
                         temp = array[m];
                         array[m] = array[i];
                         array[i] = temp;
                         break;
                    }
                    
               }
     
          }
    }
    
}
int rowReduction(double *array[],int row,int column)
{
     /**
      * Does: Reduces a square matrix m to an echelon form using elementary matrix
      * operations.
      * Recieves:
      * @array: an array of pointers, where each pointer points to a row of an augmented matrix.
      * @row: the number of rows in the augumented matrix.
      * @column: the number of columns in the coefficient matrix.
      * Returns: 1 if successful,otherwise it returns a zero
     */
    double *temp;
    for (int m = 0; m < (row - 1); m++)
     {
           if (array[m][m] == 0 )
          {
               /*scan through column n, and make a swap if there is a non-zero row*/
               for (int i = m; i < row; i++)
               {
                    if (array[i][m] > 0)
                    {
                         /*Make the swap*/
                         temp = array[m];
                         array[m] = array[i];
                         array[i] = temp;
                         break;
                    }
                    
               }
     
          }
               /*search for the first nonzero entry in the row*/
                    int n = 0;
                    for (n = 0; n < column; n++)
                    {
                         if (array[m][n] != 0)
                         {
                             
                              break;
                         }
                         
                    }
                    if (n >= column)
                    {
                         
                         continue;
                    }
                    
               /**scan the column again, starting from the next row */
               /*if their is a non-zero entry,multiply the row  and add it to that row*/
               for (int i = (m+1); i < row; i++)
               {
                    
                    double multiplier = ((array[i][n])/(array[m][n])*(-1.0));
                    
                    if (multiplier != 0.0)
                    {
                         for (int j = m; j < (column + 1); j++)
                         {

                              *(array[i] + j) = (*(array[m] + j))*multiplier + *(array[i] + j);
                              

                         }
                   
                    }
                    
               } 
  
          
     }

     return (1);

}
void gaussJordan(double *echelon[],int row,int column)
{
     /**
      * Does: Reduces an echelon matrix into a Jordan matrix
      * Recieves:
      * @echelon: an array of pointers to rows of the augmented matrix
      * @row: the number of rows in the augumented matrix.
      * @column: the number of columns in the coefficient matrix.
     */
     for (int i = 0; i < row; i++)
     {
          /*search for the first nonzero entry in the row*/
           int n = 0;
          for (n = 0; n < column; n++)
          {
               if (echelon[i][n] != 0)
               {
                    break;
               }
               
          }
          if (n >= column)
          {
                         
                continue;
          }
          double divisor = *(echelon[i] + n);
          
          for (int j = 0; j < (column + 1); j++)
          {
              *(echelon[i] + j) =  *(echelon[i] + j)/divisor;
          }
          
     }
     for (int m = (row - 1); m > 0; m--)
     {
          
          /**scan the column again, starting from the next row */
          /*if their is a non-zero entry,multiply the row  and add it to that row*/
          /**
           * First search for the first non-zero element in this row.
          */
         int n = 0;
          for (n = 0; n < column; n++)
          {
               if (echelon[m][n] != 0)
               {
                    break;
               }
               
          }
           if (n >= column)
          {
                         
               continue;
          }
          for (int i = (m - 1); i >= 0; i--)
          {
               double multiplier = ((echelon[i][n])/(echelon[m][n])*(-1.0));
               if (multiplier != 0.0)
               {
                   
                    for (int j = n; j < (column + 1 ); j++)
                    {

                         *(echelon[i] + j) = (*(echelon[m] + j))*multiplier + *(echelon[i] + j);
                         

                    }
               
               }
               
          } 
  
          
     }
     
     
}
double* particularSoln(double *reducedEchelon[], int row,int column)
{
    /**
     * Does: extracts out the particular solution from the reduced echelon form
     * Returns: a pointer to an array containing the solutions in their order.
    */
   /*Allocate space for a column x 1 matrix*/
   double *result = (double *)malloc(sizeof(double)*(column));
   for (int i = 0; i < column; i++)
   {
     /* Initialise the data to all zeros */
     *(result + i ) = 0.0;
    
   }
   
   for (int i = 0; i < row; i++)
   {
          /* scan each ros for the first non zero entry*/
          for (int j = 0; j < column; j++)
          {
               if (reducedEchelon[i][j] > 0)
               {
                   *(result + j) = *(reducedEchelon[i] + column);
               }
               
          }
          
   }
   
     return (result);
         
}
int main(void)
{
     double *array[4]; /* An array of three pointers each to a given row.*/
     double row1[] = {1,3,0,0,65};
     double row4[] = {1,2,0,0,64};
     double row2[] = {1,0,0,0,0};
     double row3[] = {1,0,0,0,0};
     array[0] = row1;
     array[1] = row2;
     array[2] = row3;
     array[3] = row4;
     
     swapRows(array,4,4);
     rowReduction(array,4,4);
     gaussJordan(array,4,4); 
     double *solns = particularSoln(array,4,4);
     printf("Solutions: [ ");
     for (int i = 0; i < 4; i++)
     {
          printf("%5.2f ",solns[i]);
     }
     printf("]\n");
     free(solns);
     return (0);
}


