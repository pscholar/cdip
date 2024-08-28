#ifndef AFFINETRANSFORMS_H
#define AFFINETRANSFORMS_H

typedef struct 
{ 
  double y; 
  double x;
  double z;
} AffinePoint;

typedef struct
{
  AffinePoint topleft;
  AffinePoint topright;
  AffinePoint botmleft;
  AffinePoint botmright;
} AffineRect;

int
affine_compute_transrect(double tr_matrix[][3], AffineRect *transrect, int inwidth, int inheight);

int
affine_compute_boundrect(AffineRect *boundrect, AffineRect *transrect);

int
affine_compute_outdim(AffineRect *boundrect, int *outwidth, int *outheight, double grid_interval);

void **
affine_compute_transform_image(const void **image, AffineRect *boundrect, double inv_matrix[][3],int inwidth, 
      int inheight, int outwidth, int outheight, int nchannels, int nb,double grid_interval);

int
affine_compute_inv_matrix3x3(double tr_matrix[][3], double inv_matrix[][3],double *det);

void **
affine_transform_image(const void **image, double tr_matrix[][3], int inwidth, int inheight,
     int nchannels, int *outwidth, int *outheight, int nb ,double grid_interval);

#endif