#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "affinetransforms.h"
#include "image_reader_writer.h"

static int
affine_map_point(int xlook, int ylook, int *xget, int *yget, int width, int height);

void **
affine_transform_image(const void **image, double tr_matrix[][3], int inwidth, int inheight,
      int nchannels, int *outwidth, int *outheight, int nb ,double grid_interval)
{
  AffineRect transrect;
  AffineRect boundrect;
  double inverse_matrix[3][3];
  double det;
  int retval;

  retval = affine_compute_inv_matrix3x3(tr_matrix, inverse_matrix, &det);
  if (retval < 0) { return NULL;}
  retval = affine_compute_transrect(tr_matrix,&transrect,inwidth, inheight);
  if (retval < 0) { return NULL;}
  retval = affine_compute_boundrect(&boundrect, &transrect);
  if (retval < 0) { return NULL;}
  retval = affine_compute_outdim(&boundrect,outwidth, outheight, grid_interval);
  if (retval < 0) { return NULL;}
  return (affine_compute_transform_image(image, &boundrect, inverse_matrix,
          inwidth,inheight, *outwidth, *outheight,nchannels,nb, grid_interval));
}

int
affine_compute_transrect(double tr_matrix[][3], AffineRect *transrect, int inwidth, int inheight)
{
  double vec_matrix[3][4] = {
    {0, inwidth - 1,0,inwidth - 1},
    {0,0,inheight - 1,inheight - 1},
    {1,1,1,1}
  };

  double vec_trans[3][4] = {0};
  double prdt_sum;

  if(!transrect){return -1;}

  for (int k = 0; k < 4; k++)
  {
    for (int i = 0; i < 3; i++)
    {
      prdt_sum = 0;
      for (int j = 0; j < 3; j++)
      {
        prdt_sum += tr_matrix[i][j] * vec_matrix[j][k];
      }
      vec_trans[i][k] = prdt_sum;
    }
  }

  for (int i = 0; i < 4; i++)
  {
    if(vec_trans[2][i] == 0.0) { return -1;}
    for (int j = 0; j < 3; j++)
    {
      vec_trans[j][i] /= vec_trans[2][i]; 
    }
  }

  transrect->topleft.x = vec_trans[0][0];
  transrect->topleft.y = vec_trans[1][0];
  transrect->topleft.z = vec_trans[2][0];
  transrect->topright.x = vec_trans[0][1];
  transrect->topright.y = vec_trans[1][1];
  transrect->topright.z = vec_trans[2][1];  
  transrect->botmleft.x = vec_trans[0][2];
  transrect->botmleft.y = vec_trans[1][2];
  transrect->botmleft.z = vec_trans[2][2];
  transrect->botmright.x = vec_trans[0][3];
  transrect->botmright.y = vec_trans[1][3];
  transrect->botmright.z = vec_trans[2][3];

  return 1;
}

int
affine_compute_boundrect(AffineRect *boundrect, AffineRect *transrect)
{
  if((!boundrect) || (!transrect)) {return -1;}
  double xmin, ymin, xmax, ymax;
  double xarr[4] = {
    transrect->topleft.x, 
    transrect->topright.x, 
    transrect->botmleft.x,
    transrect->botmright.x
  };
  double yarr[4] = {
    transrect->topleft.y,
    transrect->topright.y,
    transrect->botmleft.y,
    transrect->botmright.y
  };

  xmin = xarr[0]; xmax = xarr[0];
  ymin = yarr[0]; ymax = yarr[0];

  for (int i = 0; i < 4; i++)
  {
    xmin = xarr[i] < xmin ? xarr[i] : xmin;
    ymin = yarr[i] < ymin ? yarr[i] : ymin;
    xmax = xarr[i] > xmax ? xarr[i] : xmax;
    ymax = yarr[i] > ymax ? yarr[i] : ymax;
  }
  
  boundrect->topleft.x  = xmin;
  boundrect->topleft.y  = ymin;
  boundrect->topleft.z  = 1;
  boundrect->topright.x = xmax;
  boundrect->topright.y = ymin;
  boundrect->topright.z = 1;
  boundrect->botmleft.x = xmin;
  boundrect->botmleft.y = ymax;
  boundrect->botmleft.z = 1;
  boundrect->botmright.x = xmax;
  boundrect->botmright.y = ymax;
  boundrect->botmright.z = 1;
  // printf("\n========Bounding AffineRectangle==========\n");
  // affine_print_rect(boundrect);
  return 1;
}

int
affine_compute_outdim(AffineRect *boundrect, int *outwidth, int *outheight, double grid_interval)
{
  if((!boundrect) || (!outwidth) || (!outheight) || (!grid_interval)) { return -1; }

  *outwidth = (boundrect->topright.x - boundrect->topleft.x) / grid_interval + 1;
  *outheight = (boundrect->botmleft.y - boundrect ->topleft.y) / grid_interval + 1; 
  return 1;
}

void **
affine_compute_transform_image(const void **image, AffineRect *boundrect, double inv_matrix[][3],int inwidth, 
      int inheight, int outwidth, int outheight, int nchannels ,int nb,double grid_interval)
{
  if((!image) || (!boundrect)){ return NULL; }

  double prdt_sum;
  double vector[3], trvector[3];
  int ox, oy, py, px,isvalid;
  vector[2] = 1;

  unsigned char **img = (unsigned char **)image;
  unsigned char **outimg = NULL;
  size_t row_stride = outwidth * nchannels * nb;
  size_t num_image_bytes = outheight * row_stride;
  outimg = (unsigned char **)img_get_row_pointers(malloc(num_image_bytes),
                                outwidth,outheight,nchannels,nb);
  if(!outimg) { return NULL;}
  
  for (int y = 0; y < outheight; y++)
  {  
    vector[1] = (y * grid_interval) + boundrect->topleft.y;
    for ( int x = 0; x < outwidth; x++)
    {
      vector[0] = (x * grid_interval) + boundrect->topleft.x;    
      for (int i = 0; i < 3; i++)
      {
        prdt_sum = 0;
        for (int j = 0; j < 3; j++)
        {
          prdt_sum += inv_matrix[i][j] * vector[j];
        }
        trvector[i] = prdt_sum;
      }
      if(trvector[2] == 0)
      {
        return (NULL);
        // TODO: free up allocated image and return NULL
      }
      trvector[0] /= trvector[2];
      trvector[1] /= trvector[2];
      ox = (int)trvector[0];
      oy = (int)trvector[1];

      /*search for nearest neighbor of ox and oy*/
      px = py = 65536;
      isvalid = affine_map_point(ox,oy,&px,&py,inwidth,inheight);
      
      if (isvalid < 0)
      {
        for (int k = 0; k < nchannels; k++)
        {
          memset((outimg[y] + (x * nchannels + k) * nb),0,nb);          
        }  
      }
      else
      {
        for (int k = 0; k < nchannels; k++)
        {
          memcpy(outimg[y] + (x * nchannels + k) * nb, img[py] + (px * nchannels + k) * nb, (size_t)nb);
        }  
      }         
    }   
  }
  
  return (void **)outimg; 
}

static int
affine_map_point(int xlook, int ylook, int *xget, int *yget, int width, int height)
{
  if((xlook >= 0 && xlook < width) && (ylook >= 0 && ylook < height)) 
  {
    *xget = xlook;
    *yget = ylook;
    return 1;
  }

  int xleft, ytop, isvalidx, isvalidy;
  xleft = xlook - 1;
  ytop = ylook - 1;

  for (int i = xleft; i < xleft + 3; i++)
  {
    isvalidx = i >= 0 && i < width ? 1 : 0;
    for(int j = ytop; j < ytop + 3; j++)
    {
      isvalidy = isvalidx && j >= 0 && j < height ? 1 : 0;
      if(isvalidy) { break; }
    }
    if(isvalidy) 
    { 
      *xget = i;
      *yget = i; 
      return 1;
    }
  }
  
  return -1;

}

int
affine_compute_inv_matrix3x3(double tr_matrix[][3], double inv_matrix[][3], double *det)
{
  double matrix_adj[3][3];

  if(!det) {return -1;}

  matrix_adj[0][0] = tr_matrix[1][1] * tr_matrix[2][2] - tr_matrix[2][1] * tr_matrix[1][2];
  matrix_adj[1][0] = -(tr_matrix[1][0] * tr_matrix[2][2] - tr_matrix[2][0] * tr_matrix[1][2]);
  matrix_adj[2][0] = tr_matrix[1][0] * tr_matrix[2][1] - tr_matrix[2][0] * tr_matrix[1][1];

  *det = tr_matrix[0][0] * matrix_adj[0][0] + tr_matrix[0][1] * matrix_adj[1][0] 
        +  tr_matrix[0][2] * matrix_adj[2][0];
  if(!(*det)) { return -1; }

  matrix_adj[0][1] = -(tr_matrix[0][1] * tr_matrix[2][2] - tr_matrix[2][1] * tr_matrix[0][2]);
  matrix_adj[1][1] = tr_matrix[0][0] * tr_matrix[2][2] - tr_matrix[2][0] * tr_matrix[0][2];
  matrix_adj[2][1] = -(tr_matrix[0][0] * tr_matrix[2][1] - tr_matrix[2][0] * tr_matrix[0][1]);

  matrix_adj[0][2] = tr_matrix[0][1] * tr_matrix[1][2] - tr_matrix[1][1] * tr_matrix[0][2];
  matrix_adj[1][2] = -(tr_matrix[0][0] * tr_matrix[1][2] - tr_matrix[1][0] * tr_matrix[0][2]);
  matrix_adj[2][2] = tr_matrix[0][0] * tr_matrix[1][1] - tr_matrix[1][0] * tr_matrix[0][1];
 
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      inv_matrix[j][i] = matrix_adj[j][i] / (*det);
    }    
  }
  return 1;
}



