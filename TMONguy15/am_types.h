/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* am_types.h                  S. Paine rev. 2019 February 27
*
* Defined types for am.
************************************************************/

#ifndef AM_AM_TYPES_H
#define AM_AM_TYPES_H

#include <stdio.h>

typedef struct model_t model_t;
typedef struct layer_t layer_t;
typedef struct column_t column_t;
typedef struct abscoeff_t abscoeff_t;
typedef struct cat_entry_t cat_entry_t;
typedef struct line_coupling_table_entry_t line_coupling_table_entry_t;
typedef struct simplex_t simplex_t;
typedef struct fit_data_t fit_data_t;

extern const model_t    MODEL_INIT;
extern const layer_t    LAYER_INIT;
extern const column_t   COLUMN_INIT;
extern const abscoeff_t ABSCOEFF_INIT;
extern const simplex_t  SIMPLEX_INIT;
extern const fit_data_t FIT_DATA_INIT;

/*
 * gridsize_t is the integer type used for indexing elements of arrays
 * associated with the frequency grid.  This includes spectra and
 * arrays used to compute their numeric transforms.
 *
 * gridsize_t must be a signed integral type for compatibility with
 * OpenMP implementations before version 3.0, which require a signed
 * loop index on parallel loops.  A 32-bit signed integer is
 * sufficient for most practical work on 64-bit as well as 32-bit
 * platforms, and will often result in slightly faster code.
 */
#define GRIDSIZE_T_MAX INT_MAX
typedef int gridsize_t;

/*
 * gridrange_t is a type used to hold a range [min, max] of grid
 * indices.
 */
typedef struct gridrange_t {
    gridsize_t min;
    gridsize_t max;
} gridrange_t;

/*
 * The model_t data structure is the top-level data structure for the
 * am model.  It contains an array of layer data structures, other
 * model configuration data, and pointers to arrays containing model
 * output spectra.
 *
 * When repeated model computations are being carried out during a
 * fit or when computing numerical Jacobians, a second model_t
 * structure (called lmodel, for "last model", throughout the program)
 * contains a copy of all the model state variables and layer data
 * from the most recently-computed model.  lmodel is used to work out
 * which parts of the model need to be recomputed after a change,
 * avoiding redundant computations.  The spectral arrays in lmodel are
 * not needed, so the corresponding pointers are set to NULL.
 */
struct model_t {
    double fmin;            /* minimum grid frequency   [GHz]               */
    double fmax;            /* maximum grid frequency   [GHz]               */
    double df;              /* frequency grid interval  [GHz]               */
    double fout_min;        /* minimum output frequency [GHz]               */
    double fout_max;        /* maximum output frequency [GHz]               */
    double flo;             /* LO frequency for IF spectra [GHz]            */
    double fif_0;           /* lowest grid frequency for IF spectra [GHz]   */
    double fif_delta;       /* offset between model and IF grid     [GHz]   */
    double fif_min;         /* minimum IF frequency [GHz]                   */
    double fif_max;         /* maximum IF frequency [GHz]                   */
    double ils_fwhm;        /* FWHM of the instrumental line shape   [GHz]  */
    double ils_fif;         /* IF frequency for heterodyne ILS modes [GHz]  */
    double dsb_utol_ratio;  /* DSB upper/lower sideband ratio               */
    double rx_gain_factor;  /* receiver gain correction factor              */
    double RH_offset_exp;   /* RH offset, stored as exp(RH_offset)          */
    double RH_scale;        /* RH scale factor                              */
    double g;               /* gravitational constant at R0 [m s^-2]        */
    double dg_dz;           /* vertical gravity gradient at R0 [s^-2]       */
    double R0;              /* planetary surface radius [cm]                */
    double z0;              /* height of base layer relative to R0 [cm]     */
    double p;               /* refract. invariant n*r*sin(za) or n*sin(za)  */
    double Pobs;            /* Pressure at observing level [mbar]           */
    double zobs;            /* User-defined height at observing level [cm]  */
    double Psource;         /* Pressure at source level [mbar]              */
    double zsource;         /* User-defined height at source level [cm]     */
    double Ptan;            /* Pressure at tangent level [mbar]             */
    double ztan;            /* User-defined height at tangent level [cm]    */
    double tol;             /* fast linesum tolerance                       */
    double selfbroad_vmr_tol; /* vmr deviation tol in selfbroad computation */
    double T0;              /* background temperature      [K]              */
    double Trx;             /* receiver noise temperature  [K]              */
    double Tref;            /* reference temperature for Y [K]              */
    double kcache_Tmin;     /* Cacheable temperature range [K], and         */
    double kcache_Tmax;     /*  interpolation interval [K], for             */
    double kcache_dT;       /*  absorption coefficients.                    */
    double sec_za;          /* secant of the zenith angle at obs. point     */
    double za;              /* zenith angle at obs. point                   */
    double am_runtime;      /* atmospheric model run time                   */
    double od_runtime;      /* optical depth computations run time          */
    double rt_runtime;      /* radiative transfer run time                  */
    double spec_runtime;    /* derived spectra run time                     */
    double runtime;         /* total computation run time                   */
    double *f;              /* frequency grid           [GHz]               */
    double *f2;             /* squared frequency grid [GHz^2]               */
    double *fif;            /* IF frequency grid        [GHz]               */
    double *ils;            /* instrumental line shape (HT, bit reversed)   */
    double *ilsworkspace;   /* workspace for convolutions with ils          */
    double *tau;            /* spectral opacity (optical depth)             */
    double *tx;             /* spectral transmittance                       */
    double *I0;             /* initial background radiance spectrum         */
    double *I;              /* spectral radiance [watt / (cm^2 * GHz * sr)] */
    double *I_ref;          /* spectral radiance at Tref                    */
    double *I_diff;         /* spectral radiance difference I - Iref        */
    double *Tb;             /* spectral Planck brightness temperature   [K] */
    double *Trj;            /* spectral Rayleigh-Jeans brightness temp. [K] */
    double *Tsys;           /* (Trj + Trx) * rx_gain_factor                 */
    double *Y;              /* [(Trj + Trx) / (Tref + Trx)] * rx_gain_factor*/
    double *L;              /* computed delay spectrum                      */
    layer_t **layer;        /* array of pointers to layer data structures   */
    gridsize_t ngrid;       /* number of frequency grid points              */
    gridsize_t nif;         /* number of IF spectrum grid points            */
    gridsize_t npad;        /* padded length of arrays for ILS convolution  */
    gridsize_t nLpad;       /* padded length of arrays for delay computation*/
    gridrange_t ilsb;       /* model grid range corresponding to lsb        */
    gridrange_t iusb;       /* model grid range corresponding to usb        */
    gridrange_t isub[2];    /* model subgrid ranges for output spectra      */
    int nlayers;            /* number of layers                             */
    int path_begin;         /* layer index for start of downward path seg.  */
    int path_mid;           /* end of downward & start of upward path seg.  */
    int path_end;           /* layer index for end of upward path segment   */
    int nkcache;            /* max num cached abs. coeff. arrays per col.   */
    int fmin_unitnum;       /* user unit index number for fmin              */
    int fmax_unitnum;       /* user unit index number for fmax              */
    int df_unitnum;         /* user unit index number for df                */
    int fout_min_unitnum;   /* user unit index number for fout_min          */
    int fout_max_unitnum;   /* user unit index number for fout_max          */
    int flo_unitnum;        /* user unit index number for flo               */
    int fif_min_unitnum;    /* user unit index number for fif_min           */
    int fif_max_unitnum;    /* user unit index number for fif_max           */
    int g_unitnum;          /* user unit index number for g                 */
    int dg_dz_unitnum;      /* user unit index number for dg_dz             */
    int R0_unitnum;         /* user unit index number for R0                */
    int z0_unitnum;         /* user unit index number for z0                */
    int Pobs_unitnum;       /* user unit index number for Pobs              */
    int zobs_unitnum;       /* user unit index number for zobs              */
    int Psource_unitnum;    /* user unit index number for Psource           */
    int zsource_unitnum;    /* user unit index number for zsource           */
    int Ptan_unitnum;       /* user unit index number for Ptan              */
    int ztan_unitnum;       /* user unit index number for ztan              */
    int geometry;           /* geometry mode and status bits                */
    int ifmode;             /* IF spectrum computation mode                 */
    int ils_fwhm_unitnum;   /* user unit index number for ils_fwhm          */
    int ils_fif_unitnum;    /* user unit index number for ils_fif           */
    int ilsmode;            /* ILS mode (normal, dsb, ssb)                  */
    int ils_typenum;        /* instrumental line shape type number          */
    int RH_adj_flag;        /* set if RH_offset or RH_gain have been set    */
    int T0_unitnum;         /* user unit index number for T0                */
    int Tref_unitnum;       /* user unit index number for Tref              */
    int Trx_unitnum;        /* user unit index number for Trx               */
    int PTmode;             /* P,T mode, see comment below                  */
    int za_unitnum;         /* user unit index number for zenith angle      */
    int log_runtimes;       /* flag to log run times and report in output   */
    int headers;            /* flag to include column headers in output     */
};

/*
 * Layers
 *
 * The layer_t data structure holds configuration data for each layer
 * in the model, including the pressure and temperature variables
 * discussed above, the lineshape types to be used on the layer, and
 * an array of column structures describing each absorbing species
 * contained within the layer.  It also contains a pointer to arrays
 * containing the layer Planck spectrum, opacity spectrum, and
 * transmittance spectrum.
 */

struct layer_t {
    double P;           /* Pressure [mbar]                                  */
    double T;           /* Temperature [K]                                  */
    double dP;          /* differential base pressure [mbar]                */
    double dphi;        /* path central angle for spherical geometry        */
    double Pbase;       /* base pressure [mbar]                             */
    double Tbase;       /* base temperature [K]                             */
    double h;           /* layer thickness [cm]                             */
    double m;           /* airmass                                          */
    double Mair;        /* average molar mass on layer                      */
    double M0;          /* default molar mass of air on this layer          */
    double gbase;       /* gravitational constant at zbase                  */
    double zbase;       /* base level height in hydrostatic models          */
    double za_base;     /* zenith angle at layer base level                 */
    double *B;          /* Planck function at T or Tbase (dep. on PTmode)   */
    double *tau;        /* line-of-sight optical depth for the layer        */
    double *tx;         /* line-of-sight transmittance exp(-tau)            */
    column_t **column;  /* array of ptrs to columns within this layer       */
    int *lineshape;     /* lineshapes, indexed by k type                    */
    int *strict_selfbroad;  /* self-broadening flags, indexed by k type     */
    int *Mair_flag;     /* Mair flags, indexed by column type               */
    double runtime;     /* run time for this layer                          */
    int tagnum;         /* index for layer tag string                       */
    int type;           /* layer type - see enum above                      */
    int ncols;          /* number of columns                                */
    int P_unitnum;      /* user unit index number for pressures             */
    int T_unitnum;      /* user unit index number for temperatures          */
    int h_unitnum;      /* user unit index number for layer thickness       */
    int updateflag;     /* flag indicating arrays need to be recomputed     */
    int vmr_stat;       /* keeps track of mixing ratio computation and use  */
};

/*
 * Columns
 *
 * Layers in am contain mixtures of absorbing species organized into
 * spatially-coincident columns.  In the simplest case, a column
 * represents the column density of a particular molecular species
 * within a layer.  However, the notion of a column generalizes to
 * include composite column types such as "dry_air" that involve
 * mixtures of different species with fixed partial mixing ratios, and
 * parametric column types such as "atten" (modeling gray attenuation)
 * that are phenomenological absorption models with N serving as a
 * scaling parameter.
 *
 * Computation of column densities and mixing ratios in am is a bit
 * subtle, since there are several ways they can be defined.  In the
 * model configuration file, the column density is either defined
 * explicitly by specifying N directly, or implicitly via a mixing
 * ratio (stored as xvmr = vmr/(1 - vmr) ).  For each column type,
 * there is also a global scaling factor Nscale which multiplies N for
 * all columns of a given type.  By default, Nscale = 1.0, but it can
 * be altered in the model configuration file, or used as a fit
 * variable, to implement scaling of a density profile for a column
 * type.  N_scaled is the column density subsequent to this scaling,
 * and vmr_scaled is the recomputed mixing ratio based on the scaled
 * column density.
 *
 * unres_lines keeps a count of unresolved lines, if any, encountered
 * during computation of all line-by-line absorption coefficients
 * associated with this column.
 */
struct column_t {
    double N;           /* column density [cm^-2], or numerical parameter   */
    double N_scaled;    /* N * scale factor for col type and layer label    */
    double xvmr;        /* xvmr = vmr / (1 - vmr), maps [0, 1] to [0, inf]  */
    double vmr_scaled;  /* volume mixing ratio derived from N_scaled        */
    double RH;          /* relative humidity as a percentage (100 = 100%)   */
    double RH_adj;      /* relative humidity adjusted for offset and scale  */
    double *ztau;       /* zenith (normal to layer boundary) optical depth  */
    abscoeff_t **abscoeff; /* array of pointers to absorp. coeff. structs   */
    double runtime;     /* run time for this column                         */
    int col_typenum;    /* numerical index of column type                   */
    int n_abscoeffs;    /* number of associated absorption coefficients     */
    int N_mode;         /* column density computation mode.  See column.h.  */
    int N_unitnum;      /* user unit index number for N                     */
    int vmr_stat;       /* keeps track of mixing ratio computation and use  */
    int unres_lines;    /* total unresolved lines for all abscoeffs         */
};

/*
 * Absorption Coefficients
 *
 * For each absorption coefficient associated with a column, an
 * abscoeff_t structure contains the index of the absorption
 * coefficient type, k_typenum, and the most recently computed k
 * array.  Absorption coefficients may be molecular [cm^2], or binary
 * [cm^5].  For molecular absorption coefficients, typically
 * associated with line-by-line absorption, opacity is tau = k * N,
 * where N is a column density [cm^-2].  For binary coefficients,
 * typically associated with collision-induced spectra and continua,
 * opacity is tau = k * N * n, where n is the volume density [cm^-3]
 * of either the same gas associated with N, or of a foreign
 * perturber.
 *
 * To save time during fits, k arrays which vary with T only can be
 * computed on a fixed T grid, as needed, and cached to memory.  k
 * arrays at intermediate T values are then interpolated between grid
 * values of T.  kcache is the array of cached k arrays.

 * vmr_selfbroad relates to the self-broadening contribution to
 * collisional linewidths in line-by-line absorption coefficient
 * computations.  Line-by-line absorption coefficients are functions
 * of mixing ratio because the self- and air-broadening coefficients
 * are different.  However, in typical applications, species like O2
 * with signicant self-broadening have essentially constant mixing
 * ratio.   Species like O3 and H2O have highly variable mixing
 * ratios, but they are generally small enough that self- broadening
 * can be ignored.  Consequently, as an optimization to minimize
 * recomputation of negligibly differing absorption coefficients, am
 * assumes default dry-air mixing ratios when computing collisional
 * linewidths, unless the strict_self_broadening flag has set on a
 * given layer for the a given absorption coefficient type.  In the
 * latter case, vmr_selfbroad is set to vmr_scaled.
 */
struct abscoeff_t {
    double vmr_selfbroad;   /* vmr for computing self-broadening            */
    double *k;              /* absorption coefficient ([cm^2] or [cm^5])    */
    double **kcache;        /* memory-cached absorption coefficients        */
    int k_typenum;          /* index num. of absorption coefficient type    */
    int unres_lines;        /* count of unresolved in-band lines            */
};

/*
 * Spectral line catalogs in am consist of static arrays, derived
 * mainly from the HITRAN database, of the following data structure.
 */
struct cat_entry_t {
    double freq;     /* line frequency                     [GHz]    */
    double S;        /* line intensity              [cm^2 * GHz]    */
    double Elo;      /* lower-state energy                   [K]    */
    float gam_air;   /* air-broadening coefficient  [GHz / mbar]    */
    float gam_self;  /* self-broadening coefficient [GHz / mbar]    */
    float nair;      /* T dependence exponent of gam_air            */
    float delta_air; /* air-induced pressure shift  [GHz / mbar]    */
    int iso;         /* HITRAN isotopologue number                  */
};

/*
 * For species with line mixing implemented (presently o2 only), the
 * catalog table is split into a table of uncoupled lines, and a table
 * of coupled lines.  Each row of the coupled line table has a corresponding
 * entry in a table of coupling parameters.
 *
 * Currently only air-induced coupling is implemented.  The parametrization
 * is that of
 *   D. S. Makarov, et al. 2011, "60-GHz oxygen band: Precise experimental
 *   profiles and extended absorption modeling in a wide temperature range."
 *   JQSRT 112:1420,
 * which follows
 *   E. W. Smith 1981, "Absorption and dispersion in the O2 microwave spectrum
 *   at atmospheric pressures."  JQSRT 74:6658.
 */
struct line_coupling_table_entry_t {
    float y0;   /* first-order mixing coefficient                       */
    float y1;   /*   temperature dependence linear parameter            */
    float yx;   /*   temperature dependence exponent                    */
    float g0;   /* second-order mixing coefficient on intensity         */ 
    float g1;   /*   temperature dependence linear parameter            */
    float gx;   /*   temperature dependence exponent                    */
    float d0;   /* second-order mixing coefficient on line frequency    */
    float d1;   /*   temperature dependence linear parameter            */
    float dx;   /*   temperature dependence exponent                    */
    float Tref; /* reference temperature for all entries                */
};


/*
 * A subset of model parameters can be varied to fit spectra using the
 * downhill-simplex method.  The following structure describes the
 * simplex and associated fit-related data.  The coordinates of
 * the simplex space are the user-defined fit variables, implemented
 * as aliased pointers.
 *
 * Jacobians can also be computed over a user-defined set of model
 * variables, and the same simplex structure is used to define the
 * variables of differentiation.
 */
struct simplex_t {
    double E1;          /* value of fit estimator at p1                     */
    double E2;          /* value of fit estimator at p2                     */
    double tol;         /* convergence tolerance for fits                   */
    char **name;        /* text names of variable parameters                */
    double *init;       /* initial values for variable parameters           */
    double *scale;      /* characteristic scale for each variable parameter */
    int *unitnum;       /* user unit numbers for variable parameters        */
    int *mapping;       /* mapping type for each variable parameter         */
    double **varptr;    /* array of pointers to variable parameters         */
    double **vertex;    /* array of vertex vectors                          */
    double *E;          /* value of fit estimator at each vertex            */
    double *pbar;       /* centroid of face opposite highest vertex         */
    double *p1;         /* trial vertex                                     */
    double *p2;         /* trial vertex                                     */
    double *pc;         /* lowest vertex at last convergence                */
    unsigned int ihigh; /* index of highest vertex                          */
    unsigned int ilow;  /* index of lowest vertex                           */
    unsigned int n;     /* simplex dimensionality                           */
    int iter;           /* maximum iteration count for fits                 */
    int reinit;         /* flag to reset initial values for new fit         */
    int restart_count;  /* counts restarts after convergence                */
    int logarithmic;    /* flag indicating log mapping of coordinates       */
};

/*
 * For model fits, the following structure contains the list of files
 * containing data to be fit, arrays for the current fit data, and
 * miscellaneous fit-related variables.
 *
 * The initializer for this structure is defined in fit.c.
 */
struct fit_data_t {
    char **filename;    /* array of file names                              */
    FILE *fp;           /* pointer for currently open fit data file         */
    double *f;          /* fit data frequency                               */
    double *s;          /* fit data spectrum                                */
    double *s_mod;      /* model data interpolated onto fit data grid       */
    double *res;        /* fit residuals                                    */
    double *res_est;    /* fit residuals estimated by recursive filter      */
    double *b;          /* fit data channel bandwidth                       */
    double *w;          /* fit data weight factor (1 / sigma)               */
    double mean_var;    /* mean variance or reduced chi_squared             */
    double mean_var_tracked; /* mean variance after residual tracking       */
    double res_track_gain;   /* gain for residual tracking                  */
    double runtime;     /* fit run time                                     */
    gridsize_t npts;    /* number of data points                            */
    gridsize_t nalloc;  /* currently allocated size of fit data arrays      */
    int data_type;      /* type of data being fit                           */
    int estimator_type; /* type of fit estimator being used                 */
    int f_unitnum;      /* unit index number for frequency data             */
    int s_unitnum;      /* unit index number for spectrum data              */
    int f_col;          /* input column number for frequency                */
    int s_col;          /* input column number for spectrum                 */
    int b_col;          /* input column number for channel bandwidth        */
    int w_col;          /* input column number weight factor                */
    int max_iter;       /* maximum number of downhill simplex iterations    */
    int max_restarts;   /* max number of restarts after convergence         */
    int nfiles;         /* number of files to be fit                        */
    int open_filenum;   /* array index of current fit data file             */
    int blocknum;       /* current data block within current file           */
    int output_mode;    /* output mode control, see bits enumerated above   */
};

#endif /* AM_AM_TYPES_H */
