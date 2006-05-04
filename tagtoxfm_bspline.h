/* ----------------------------- MNI Header -----------------------------------
@NAME       : tagtoxfm_bspline.h
@DESCRIPTION: Header file for tagtoxfm_bspline.cxx
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 7, 2006  John G. Sled
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define DEFAULT_REAL_RANGE 3.0

int inverse = false;
int clobber = false;
int with_linear = true;
Trans_type transform_type = TRANS_LSQ6;

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
   {NULL, ARGV_HELP, NULL, NULL,
       "Linear Transformation type. Default = -lsq6."},
   {"-lsq6", ARGV_CONSTANT, (char *) TRANS_LSQ6, (char *) &transform_type,
       "6 parameter (scale=1.0) least-squares linear transformation."},
   {"-lsq7", ARGV_CONSTANT, (char *) TRANS_LSQ7, (char *) &transform_type,
       "7 parameter (one scale) least-squares linear transformation."},
   {"-lsq9", ARGV_CONSTANT, (char *) TRANS_LSQ9, (char *) &transform_type,
       "9 parameter least-squares linear transformation."},
   {"-lsq10", ARGV_CONSTANT, (char *) TRANS_LSQ10, (char *) &transform_type,
       "10 parameter least-squares linear transformation."},
   {"-lsq12", ARGV_CONSTANT, (char *) TRANS_LSQ12, (char *) &transform_type,
       "12 parameter least-squares linear transformation."},
   {"-with_linear", ARGV_CONSTANT, (char *) TRUE, (char *) &with_linear,
      "Estimate linear transform before computing non-linear component of transform (default)."},
   {"-without_linear", ARGV_CONSTANT, (char *) FALSE, (char *) &with_linear,
   "Do not estimate linear transform."},
   {NULL, ARGV_HELP, NULL, NULL,
       "Non-linear options."},
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
					 int with_linear_flag,
					 Trans_type transform_type);

int create_grid_volume_from_example(char *likefile, Volume &grid_volume); 
void determine_domain(int n_tag_points, Real **tags_volume1,
		      Volume grid_volume, DblMat &domain);
