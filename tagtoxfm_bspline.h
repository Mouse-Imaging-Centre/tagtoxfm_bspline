/* ----------------------------- MNI Header -----------------------------------
@NAME       : tagtoxfm_bspline.h
@DESCRIPTION: Header file for tagtoxfm_bspline.cxx
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 7, 2006  John G. Sled
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define MIN_POINTS 4
#define DEFAULT_REAL_RANGE 3.0

int inverse = false;
int clobber = false;
int rigid = true;
double lambda = 0.01;      // scale invariant smoothing parameter
double distance = 10;     // distance between basis functions
double real_range[2] = { -DEFAULT_REAL_RANGE, DEFAULT_REAL_RANGE };

/* Argument table */
ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, NULL, NULL,
       "---Transformation maps volume two to volume one---"},
   {"-inverse", ARGV_CONSTANT, (char *) TRUE, (char *) &inverse,
       "Swap tags, then compute transform (default=FALSE)."},
   {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
       "Overwrite any existing xfm file."},
   {"-noclobber", ARGV_CONSTANT, (char *) FALSE, (char *) &clobber,
       "Do not overwrite any existing xfm file."},
   {"-with_rigid", ARGV_CONSTANT, (char *) TRUE, (char *) &rigid,
      "Estimate rigid body transform before computing non-linear component of transform (default)."},
   {"-without_rigid", ARGV_CONSTANT, (char *) FALSE, (char *) &rigid,
   "Do not estimate rigid body transform."},
   {"-lambda", ARGV_FLOAT, (char *) 1, (char *) &lambda, 
   "Scale invariant smoothing parameter.  The default is chosen\n to provide minimal"
   " smoothing except in regions of missing data."},
   {"-distance", ARGV_FLOAT, (char *) 1, (char *) &distance, 
   "Distance between basis functions in mm.  This parameter\n determines the overall"
   " smoothness."},
   {"-range", ARGV_FLOAT, (char *) 2, (char *) &real_range, 
    "range of values for grid transform [default: -3 3]"},
   {NULL, ARGV_END, NULL, NULL, NULL}
};


int compute_tbspline_transform_from_tags(int n_tag_points, 
					 Real **tags_volume1, 
					 Real **tags_volume2, 
					 General_transform *transform,
					 Volume grid_volume, 
					 Real distance,
					 Real lambda,
					 int rigid_flag);

int create_grid_volume_from_example(char *likefile, Volume &grid_volume); 
void determine_domain(int n_tag_points, Real **tags_volume1,
		      Volume grid_volume, DblMat &domain);
