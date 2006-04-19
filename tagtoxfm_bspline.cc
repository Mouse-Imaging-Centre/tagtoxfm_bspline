
/* tagtoxfm_bspline is a program to read a set of tag points and 
   create an approximate grid transform using tensor cubic B splines 

   Created:  April 7, 2006   John G. Sled 
   using code from tagtoxfm by Peter Neelin */

#include "TBSpline.h"  // Note: this produces odd compiler errors if 
                       //  included after the other header files

extern "C" {
#include <stdio.h>
#include <string.h>
#include <ParseArgv.h>
#include <volume_io.h>
#include <bicpl/compute_xfm.h>
#define public
#include <bicpl/trans_prototypes.h>
}

#include "tagtoxfm_bspline.h"


int main(int argc, char *argv[])
{
  char *pname, *tagfile, *likefile, *xfmfile;
  int n_volumes, n_tag_points;
  Real **tags_volume1, **tags_volume2, **temp_tags;
  General_transform transform;
  Volume like_volume, grid_volume;
  char *inverse_string, comment[512];
  FILE *fp;

  pname = argv[0];
  if (ParseArgv(&argc, argv, argTable, 0) || (argc != 4)) {
    (void) fprintf(stderr, 
    "\nUsage: %s [<options>] infile.tag like_volume.mnc outfile.xfm\n\n",
                     argv[0]);
    exit(EXIT_FAILURE);
  }
  tagfile = argv[1];
  likefile = argv[2];
  xfmfile = argv[3];


  /* Read in tag file */
  if ((open_file_with_default_suffix(tagfile,
                  get_default_tag_file_suffix(),
				     READ_FILE, ASCII_FORMAT, &fp) != OK) ||
      (input_tag_points(fp, &n_volumes, &n_tag_points, 
			&tags_volume1, &tags_volume2, 
			NULL, NULL, NULL, NULL) != OK)) {
    (void) fprintf(stderr, "%s: Error reading tag file %s\n", 
		   pname, tagfile);
    exit(EXIT_FAILURE);
  }
  (void) close_file(fp);


  /* Check number of volumes */
  if (n_volumes != 2) {
    (void) fprintf(stderr, "%s: Wrong number of volumes in %s\n", 
		   pname, tagfile);
    exit(EXIT_FAILURE);
  }

  if (n_tag_points < MIN_POINTS) {
    (void) fprintf(stderr, 
		   "%s: Need at least %d points (only %d in %s)\n", 
		   pname, MIN_POINTS, n_tag_points, tagfile);
    exit(EXIT_FAILURE);
  }

  if (distance <= 0.0) {
    (void) fprintf(stderr, 
		   "%s: knot distance must be greater than 0.0\n", pname);
    exit(EXIT_FAILURE);
  }

  /* If inverting, switch order of points */
  if (inverse) {
    temp_tags = tags_volume1;
    tags_volume1 = tags_volume2;
    tags_volume2 = temp_tags;
  }


  /* create a grid transform volume using likefile as an example */
  if (create_grid_volume_from_example(likefile, grid_volume) != EXIT_SUCCESS) {
    exit(EXIT_FAILURE);
  }

  /* set output range */
  set_volume_real_range(grid_volume, real_range[0], real_range[1]);    


   /* Compute transformation */
  if(compute_tbspline_transform_from_tags(n_tag_points, tags_volume1, 
					tags_volume2, &transform,
					  grid_volume, distance, lambda,
					  rigid)  
     != EXIT_SUCCESS) {
    exit(EXIT_FAILURE);
  }
  

  if (inverse)
    inverse_string = " with -inverse";
  else
     inverse_string = "";
  (void) sprintf(comment, " Created from tag file %s\n using tensor cubic bsplines %s.",
		 tagfile, inverse_string);

  /* Save transformation */
  if (!clobber && file_exists(xfmfile)) {
     (void) fprintf(stderr, "%s: xfm file \"%s\" already exists.\n",
                    pname, xfmfile);
     exit(EXIT_FAILURE);
  }
  if (output_transform_file(xfmfile, comment, &transform) != OK) {
     (void) fprintf(stderr, "%s: Error writing xfm file %s\n", 
                    pname, xfmfile);
     exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}


/* estimate an approximate transform from the supplied tags using 
   tensor cubic B splines and use this to produce a grid transform 

   transform is an output parameter

   returns: EXIT_SUCESS or EXIT_FAILURE
*/
int compute_tbspline_transform_from_tags(int n_tag_points, 
					 Real **tags_volume1, 
					 Real **tags_volume2, 
					 General_transform *transform,
					 Volume grid_volume, 
					 Real distance,
					 Real lambda,
					 int rigid_flag) 
{
  DblMat domain(N_DIMENSIONS,2);
  Real point[N_DIMENSIONS];
  Real voxel[N_DIMENSIONS];
  float fpoint[N_DIMENSIONS];
  TBSpline *spline[N_DIMENSIONS];
  int i, j;
  int k0, k1, k2;
  Real displacement;
  General_transform linear_transform, *v_to_w_transform,
    grid_transform, v_to_w_and_linear_transform;


  if (rigid_flag) {
    /* determine an rigid transformation using the tag points */
    compute_transform_from_tags(n_tag_points, tags_volume1, tags_volume2,
				TRANS_LSQ6, &linear_transform);

  
    /* apply linear transformation to the grid volume */
    v_to_w_transform = get_voxel_to_world_transform(grid_volume);
    concat_general_transforms(v_to_w_transform, &linear_transform,
			      &v_to_w_and_linear_transform);
    set_voxel_to_world_transform(grid_volume, &v_to_w_and_linear_transform);
			    
    /* apply linear transformation to tag set 1 */
    for (j = 0; j < n_tag_points; j++) {
      general_transform_point(&linear_transform,
			      tags_volume1[j][0], tags_volume1[j][1], 
			      tags_volume1[j][2], &point[0], &point[1], 
			      &point[2]);
      for (i = 0; i < N_DIMENSIONS; i++) {
	tags_volume1[j][i] = point[i]; 
      }
    }
  }

  /* determine suitable domain for spline basis functions in voxel coords */
  determine_domain(n_tag_points, tags_volume1, grid_volume, domain);


  /* determine coordinates for sampling of volume */
  int sizes[N_DIMENSIONS+1];
  Real origin[] = { 0.0, 0.0, 0.0};
  Real separations[N_DIMENSIONS+1];
  
  get_volume_sizes(grid_volume, sizes);
  get_volume_separations(grid_volume, separations);

  /* create tensor cubic B spline basis */
  for (i = 0; i < N_DIMENSIONS; i++) {
    spline[i] = new TBSplineVolume(domain, origin, separations,
				     &sizes[1], distance, lambda, TRUE);
  }
  
  /* fit splines to tag points */
  for (j = 0; j < n_tag_points; j++) {
    convert_world_to_voxel(grid_volume, tags_volume1[j][0], 
			   tags_volume1[j][1], tags_volume1[j][2], voxel);
    for (i = 0; i < N_DIMENSIONS; i++) {
      fpoint[i] = (float) voxel[i+1]*separations[i+1];
    }
    for (i = 0; i < N_DIMENSIONS; i++) {
      spline[i]->addDataPoint(fpoint, tags_volume2[j][i]-tags_volume1[j][i]);
    }
  }
  for (i = 0; i < N_DIMENSIONS; i++) {
    if(spline[i]->fit() == FALSE) { // fit splines to the data
      fprintf(stderr, "Fatal Error: Spline fit failed.\n");
      return EXIT_FAILURE;
    }
  }

  /* evaluate spline for each voxel in grid transform volume */
  progress_struct progress;
  initialize_progress_report(&progress, FALSE, sizes[1],
			     "Evaluating splines");
  for (k0 = 0; k0 < sizes[1]; k0++) {
    for (k1 = 0; k1 < sizes[2]; k1++) {
      for (k2 = 0; k2 < sizes[3]; k2++) {
	for (i = 0; i < N_DIMENSIONS; i++) {
	  displacement = (*(spline[i]))(k0*separations[1], 
					k1*separations[2],
					k2*separations[3]);
	  set_volume_real_value(grid_volume, i, k0, k1, k2, 0, displacement);
	}
      }
    }
    update_progress_report(&progress, (int) k0);
  }
  terminate_progress_report(&progress);

  /* create grid transform using grid transform volume */
  create_grid_transform(&grid_transform, grid_volume);

  if (rigid_flag){
    /* concatenate linear transform with grid_transform */
    concat_general_transforms(transform, &linear_transform,
			      &grid_transform);
  }
  else
  {
    *transform = grid_transform;  /* this assumes grid_transform won't be
				    explicitly deleted */
  }  
  
  
  return EXIT_SUCCESS;
}


/* determine a rectangular domain that includes tag set 1 and 
   all voxels of grid_volume using voxel coordinates

   return: domain (2 by N matrix) as an output variable */
void determine_domain(int n_tag_points, Real **tags_volume1,
		      Volume grid_volume, DblMat &domain)
{
  int sizes[N_DIMENSIONS+1];
  int i, j;
  Real voxel[N_DIMENSIONS+1];
  Real separations[N_DIMENSIONS+1];

  get_volume_sizes(grid_volume, sizes);
  get_volume_separations(grid_volume, separations);

  /* set domain to that of grid volume */
  for(i = 0; i < N_DIMENSIONS; i++)
    {
      if(separations[i+1] > 0) {
        domain(i,0) = -0.5*separations[i+1];
        domain(i,1) = (sizes[i+1]-0.5)*separations[i+1];
      }
      else {
        domain(i,1) = -0.5*separations[i+1];
        domain(i,0) = (sizes[i+1]-0.5)*separations[i+1];
      }
    }


  /* determine domain of coordinates in tag sets 1 */
  for (j = 0; j < n_tag_points; j++) {
    convert_world_to_voxel(grid_volume, tags_volume1[j][0], 
			   tags_volume1[j][1], tags_volume1[j][2],
			   voxel);
    for (i = 0; i < N_DIMENSIONS; i++) {
      if (domain(i, 0) > voxel[i+1]*separations[i+1]) {
	domain(i, 0) = voxel[i+1]*separations[i+1];    /* set new minimum */
      }
      if (domain(i, 1) < voxel[i+1]*separations[i+1]) {
	domain(i, 1) = voxel[i+1]*separations[i+1];    /* set new maximum */
      }
    }
  }
}


int create_grid_volume_from_example(char *likefile, Volume &grid_volume)
{
  int sizes[N_DIMENSIONS+1];
  General_transform *transform, *new_transform;
  Volume like_volume;
  Real cosine[N_DIMENSIONS];
  Real separations[N_DIMENSIONS+1];
  int i;
  Real v_origin[4] = {0.0, 0.0, 0.0, 0.0};
  Real w_origin[3];

  char *standard_dimension_order[3] = {MIzspace, MIyspace, MIxspace};
  char *grid_dimension_order[4] = 
    { MIvector_dimension, MIzspace, MIyspace, MIxspace };

  if (input_volume(likefile, N_DIMENSIONS, standard_dimension_order,
		   NC_UNSPECIFIED, /* data type */ FALSE,
		   NC_UNSPECIFIED, NC_UNSPECIFIED, /* min, max */
		   TRUE, &like_volume, (minc_input_options *) NULL) != OK) {
    (void) fprintf(stderr, "Error reading file %s\n", likefile);
    return EXIT_FAILURE;
  }
   
  grid_volume = create_volume(N_DIMENSIONS+1, 
			      grid_dimension_order, NC_SHORT, true, 0, 0);

  get_volume_sizes(like_volume, &sizes[1]);
  sizes[0] = 3;
  set_volume_sizes(grid_volume, sizes);

  /*  transform = get_voxel_to_world_transform(like_volume);
  copy_general_transform(transform, new_transform);
  set_voxel_to_world_transform(grid_volume, new_transform); */

  for (i = 0; i < N_DIMENSIONS; i++) {
    get_volume_direction_cosine(like_volume, i, cosine);
    set_volume_direction_cosine(grid_volume, i+1, cosine);
  }
  get_volume_separations(like_volume, &separations[1]);
  separations[0] = 1;
  set_volume_separations(grid_volume, separations);
  
  convert_voxel_to_world(like_volume, v_origin, &w_origin[0], &w_origin[1], 
			 &w_origin[2]);
  set_volume_translation(grid_volume, v_origin, w_origin);

  delete_volume(like_volume);
  alloc_volume_data(grid_volume);

  return EXIT_SUCCESS;
}


