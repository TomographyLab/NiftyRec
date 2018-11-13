/*
  Teem: Tools to process and visualize scientific data and images              
  Copyright (C) 2011, 2010, 2009  University of Chicago
  Copyright (C) 2008, 2007, 2006, 2005  Gordon Kindlmann
  Copyright (C) 2004, 2003, 2002, 2001, 2000, 1999, 1998  University of Utah

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
  The terms of redistributing and/or modifying this software also
  include exceptions to the LGPL that facilitate static linking.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this library; if not, write to Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "nrrd.h"

/*
** These kernels are the cardinal B-splines of different orders
** Using them with convolution assumes that the data has been pre-filtered
** so that the spline interpolates the original values.
*/

/* helper macros for doing abs() and remembering sign */
#define ABS_SGN(ax, sgn, x)                     \
  if (x < 0) {                                  \
    sgn = -1;                                   \
    ax = -x;                                    \
  } else {                                      \
    sgn = 1;                                    \
    ax = x;                                     \
  }

/* helper macro for listing the various members of the kernel */
#define BSPL_DECL(ord, deriv)                   \
  0,                                            \
    _bspl##ord##_sup,                           \
    _bspl##ord##d##deriv##_int,                 \
    _bspl##ord##d##deriv##_1f,                  \
    _bspl##ord##d##deriv##_Nf,                  \
    _bspl##ord##d##deriv##_1d,                  \
    _bspl##ord##d##deriv##_Nd

/* ============================= order *3* ============================= */

static double
_bspl3_sup(const double *parm) {
  AIR_UNUSED(parm);
  return 2.0;
}

/* ---------------------- order *3* deriv *0* -------------------------- */

static double
_bspl3d0_int(const double *parm) {
  AIR_UNUSED(parm);
  return 1.0;
}

/* t: tmp; ax: abs(x) */
#define BSPL3D0(ret, t, x)                     \
  if (x < 1) {                                 \
    ret = (4 + 3*(-2 + x)*x*x)/6;              \
  } else if (x < 2) {                          \
    t = (-2 + x);                              \
    ret = -t*t*t/6;                            \
  } else {                                     \
    ret = 0;                                   \
  }

static double
_bspl3d0_1d(double x, const double *parm) {
  double ax, tmp, r;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL3D0(r, tmp, ax);
  return r;
}

static float
_bspl3d0_1f(float x, const double *parm) {
  float ax, tmp, r;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL3D0(r, tmp, ax);
  return r;
}

static void
_bspl3d0_Nd(double *f, const double *x, size_t len, const double *parm) {
  double ax, tmp, r;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = x[i]; ax = AIR_ABS(ax);
    BSPL3D0(r, tmp, ax);
    f[i] = r;
  }
}

static void
_bspl3d0_Nf(float *f, const float *x, size_t len, const double *parm) {
  float ax, tmp, r;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = x[i]; ax = AIR_ABS(ax);
    BSPL3D0(r, tmp, ax);
    f[i] = r;
  }
}

NrrdKernel
_nrrdKernelBSpline3 = {
  "bspl3",
  BSPL_DECL(3, 0)
};
NrrdKernel *const
nrrdKernelBSpline3 = &_nrrdKernelBSpline3;

/* ---------------------- order *3* deriv *1* -------------------------- */

static double
_bspl3d1_int(const double *parm) {
  AIR_UNUSED(parm);
  return 0.0;
}

/* t: tmp; ax: abs(x) */
#define BSPL3D1(ret, t, x)                     \
  if (x < 1) {                                 \
    ret = (-4 + 3*x)*x/2;                      \
  } else if (x < 2) {                          \
    t = (-2 + x);                              \
    ret = -t*t/2;                              \
  } else {                                     \
    ret = 0;                                   \
  }

static double
_bspl3d1_1d(double x, const double *parm) {
  double ax, tmp, r;
  int sgn;
  AIR_UNUSED(parm);

  ABS_SGN(ax, sgn, x);
  BSPL3D1(r, tmp, ax);
  return sgn*r;
}

static float
_bspl3d1_1f(float x, const double *parm) {
  float ax, tmp, r;
  int sgn;
  AIR_UNUSED(parm);

  ABS_SGN(ax, sgn, x);
  BSPL3D1(r, tmp, ax);
  return sgn*r;
}

static void
_bspl3d1_Nd(double *f, const double *x, size_t len, const double *parm) {
  double ax, tmp, r;
  int sgn;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ABS_SGN(ax, sgn, x[i]);
    BSPL3D1(r, tmp, ax);
    f[i] = sgn*r;
  }
}

static void
_bspl3d1_Nf(float *f, const float *x, size_t len, const double *parm) {
  float ax, tmp, r;
  int sgn;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ABS_SGN(ax, sgn, x[i]);
    BSPL3D1(r, tmp, ax);
    f[i] = sgn*r;
  }
}

NrrdKernel
_nrrdKernelBSpline3D = {
  "bspl3d",
  BSPL_DECL(3, 1)
};
NrrdKernel *const
nrrdKernelBSpline3D = &_nrrdKernelBSpline3D;

/* ---------------------- order *3* deriv *2* -------------------------- */

static double
_bspl3d2_int(const double *parm) {
  AIR_UNUSED(parm);
  return 0.0;
}

/* NOTE: the tmp variable wasn't actually needed here, and this will
** likely be optimized out.  But the tmp argument to the macro is kept
** here (and the macro uses it to avoid a unused variable warning) to
** facilitate copy-and-paste for higher-order splines
*/
#define BSPL3D2(ret, tmp, x)                   \
  if (x < 1) {                                 \
    ret = -2 + 3*x;                            \
  } else if (x < 2) {                          \
    tmp = 2 - x;                               \
    ret = tmp;                                 \
  } else {                                     \
    ret = 0;                                   \
  }

static double
_bspl3d2_1d(double x, const double *parm) {
  double ax, tmp, r;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL3D2(r, tmp, ax);
  return r;
}

static float
_bspl3d2_1f(float x, const double *parm) {
  float ax, tmp, r;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL3D2(r, tmp, ax);
  return r;
}

static void
_bspl3d2_Nd(double *f, const double *x, size_t len, const double *parm) {
  double ax, tmp, r;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = AIR_ABS(x[i]);
    BSPL3D2(r, tmp, ax);
    f[i] = r;
  }
}

static void
_bspl3d2_Nf(float *f, const float *x, size_t len, const double *parm) {
  float ax, tmp, r;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = AIR_ABS(x[i]);
    BSPL3D2(r, tmp, ax);
    f[i] = r;
  }
}

NrrdKernel
_nrrdKernelBSpline3DD = {
  "bspl3dd",
  BSPL_DECL(3, 2)
};
NrrdKernel *const
nrrdKernelBSpline3DD = &_nrrdKernelBSpline3DD;

/* ------------- order *3* approximate numerical inverse -------------- */
/* still need to implement:
**   Unser et al B-Spline Signal Processing: Part I & II, IEEE
**   Transactions on Signal Processing, 1993, 41(2):821-833, 834--848
** but until then here's a slower way of approximating the prefiltering,
** which is still faster than doing iterative deconvolution.  These
** weights were determined by GLK with Mathematica, by inverting the
** matrix representing discrete convolution with the spline
*/

static double
_bspl3_ANI_sup(const double *parm) {
  AIR_UNUSED(parm);
  return 12.5;
}

static double
_bspl3_ANI_int(const double *parm) {
  AIR_UNUSED(parm);
  return 1.0;
}

static double
_bspl3_ANI_kvals[12] = {
  2672279.0/1542841.0,
  -(716035.0/1542841.0),
  191861.0/1542841.0,
  -(51409.0/1542841.0),
  13775.0/1542841.0,
  -(3691.0/1542841.0),
  989.0/1542841.0,
  -(265.0/1542841.0),
  71.0/1542841.0,
  -(19.0/1542841.0),
  5.0/1542841.0,
  -(1.0/1542841.0)};

#define BSPL3_ANI(ret, tmp, x)                  \
  tmp = AIR_CAST(unsigned int, x+0.5);          \
  if (tmp < 12) {                               \
    ret = _bspl3_ANI_kvals[tmp];                \
  } else {                                      \
    ret = 0.0;                                  \
  }

static double
_bspl3_ANI_1d(double x, const double *parm) {
  double ax, r; int tmp;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL3_ANI(r, tmp, ax);
  return r;
}

static float
_bspl3_ANI_1f(float x, const double *parm) {
  double ax, r; int tmp;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL3_ANI(r, tmp, ax);
  return AIR_CAST(float, r);
}

static void
_bspl3_ANI_Nd(double *f, const double *x, size_t len, const double *parm) {
  double ax, r; int tmp;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = x[i]; ax = AIR_ABS(ax);
    BSPL3_ANI(r, tmp, ax);
    f[i] = r;
  }
}

static void
_bspl3_ANI_Nf(float *f, const float *x, size_t len, const double *parm) {
  double ax, r; int tmp;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = x[i]; ax = AIR_ABS(ax);
    BSPL3_ANI(r, tmp, ax);
    f[i] = AIR_CAST(float, r);
  }
}

NrrdKernel
_nrrdKernelBSpline3ApproxInverse = {
  "bspl3ai", 0,
  _bspl3_ANI_sup, _bspl3_ANI_int,
  _bspl3_ANI_1f, _bspl3_ANI_Nf,
  _bspl3_ANI_1d, _bspl3_ANI_Nd
};
NrrdKernel *const
nrrdKernelBSpline3ApproxInverse = &_nrrdKernelBSpline3ApproxInverse;

/* ============================= order *4* ============================= */

/*
static double
_bspl4_sup(const double *parm) {
  AIR_UNUSED(parm);
  return 2.5;
}
*/

/* ============================= order *5* ============================= */

static double
_bspl5_sup(const double *parm) {
  AIR_UNUSED(parm);
  return 3.0;
}

/* ---------------------- order *5* deriv *0* -------------------------- */

static double
_bspl5d0_int(const double *parm) {
  AIR_UNUSED(parm);
  return 1.0;
}

#define BSPL5D0(ret, t, x)                                      \
  if (x < 1) {                                                  \
    t = x*x;                                                    \
    ret = (33 - 5*t*(6 + (x-3)*t))/60;                          \
  } else if (x < 2) {                                           \
    ret = (51 + 5*x*(15 + x*(-42 + x*(30 + (-9 + x)*x))))/120;  \
  } else if (x < 3) {                                           \
    t = x - 3;                                                  \
    ret = -t*t*t*t*t/120;                                       \
  } else {                                                      \
    ret = 0;                                                    \
  }

static double
_bspl5d0_1d(double x, const double *parm) {
  double ax, tmp, r;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL5D0(r, tmp, ax);
  return r;
}

static float
_bspl5d0_1f(float x, const double *parm) {
  float ax, tmp, r;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL5D0(r, tmp, ax);
  return r;
}

static void
_bspl5d0_Nd(double *f, const double *x, size_t len, const double *parm) {
  double ax, tmp, r;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = x[i]; ax = AIR_ABS(ax);
    BSPL5D0(r, tmp, ax);
    f[i] = r;
  }
}

static void
_bspl5d0_Nf(float *f, const float *x, size_t len, const double *parm) {
  float ax, tmp, r;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = x[i]; ax = AIR_ABS(ax);
    BSPL5D0(r, tmp, ax);
    f[i] = r;
  }
}

NrrdKernel
_nrrdKernelBSpline5 = {
  "bspl5",
  BSPL_DECL(5, 0)
};
NrrdKernel *const
nrrdKernelBSpline5 = &_nrrdKernelBSpline5;

/* ---------------------- order *5* deriv *1* -------------------------- */

static double
_bspl5d1_int(const double *parm) {
  AIR_UNUSED(parm);
  return 0.0;
}

#define BSPL5D1(ret, t, x)                              \
  if (x < 1) {                                          \
    t = x*x*x;                                          \
    ret = -x + t - (5*t*x)/12;                          \
  } else if (x < 2) {                                   \
    ret = (15 + x*(-84 + x*(90 + x*(-36 + 5*x))))/24;   \
  } else if (x < 3) {                                   \
    t = -3 + x;                                         \
    ret = -t*t*t*t/24;                                  \
  } else {                                              \
    ret = 0;                                            \
  }

static double
_bspl5d1_1d(double x, const double *parm) {
  double ax, tmp, r;
  int sgn;
  AIR_UNUSED(parm);

  ABS_SGN(ax, sgn, x);
  BSPL5D1(r, tmp, ax);
  return sgn*r;
}

static float
_bspl5d1_1f(float x, const double *parm) {
  float ax, tmp, r;
  int sgn;
  AIR_UNUSED(parm);

  ABS_SGN(ax, sgn, x);
  BSPL5D1(r, tmp, ax);
  return sgn*r;
}

static void
_bspl5d1_Nd(double *f, const double *x, size_t len, const double *parm) {
  double ax, tmp, r;
  int sgn;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ABS_SGN(ax, sgn, x[i]);
    BSPL5D1(r, tmp, ax);
    f[i] = sgn*r;
  }
}

static void
_bspl5d1_Nf(float *f, const float *x, size_t len, const double *parm) {
  float ax, tmp, r;
  int sgn;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ABS_SGN(ax, sgn, x[i]);
    BSPL5D1(r, tmp, ax);
    f[i] = sgn*r;
  }
}

NrrdKernel
_nrrdKernelBSpline5D = {
  "bspl5d",
  BSPL_DECL(5, 1)
};
NrrdKernel *const
nrrdKernelBSpline5D = &_nrrdKernelBSpline5D;

/* ---------------------- order *5* deriv *2* -------------------------- */

static double
_bspl5d2_int(const double *parm) {
  AIR_UNUSED(parm);
  return 0.0;
}

#define BSPL5D2(ret, t, x)                      \
  if (x < 1) {                                  \
    t = x*x;                                    \
    ret = -1 + 3*t - (5*t*x)/3;                 \
  } else if (x < 2) {                           \
    ret = (-21 + x*(45 + x*(-27 + 5*x)))/6;     \
  } else if (x < 3) {                           \
    t = -3 + x;                                 \
    ret = -t*t*t/6;                             \
  } else {                                      \
    ret = 0;                                    \
  }

static double
_bspl5d2_1d(double x, const double *parm) {
  double ax, tmp, r;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL5D2(r, tmp, ax);
  return r;
}

static float
_bspl5d2_1f(float x, const double *parm) {
  float ax, tmp, r;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL5D2(r, tmp, ax);
  return r;
}

static void
_bspl5d2_Nd(double *f, const double *x, size_t len, const double *parm) {
  double ax, tmp, r;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = AIR_ABS(x[i]);
    BSPL5D2(r, tmp, ax);
    f[i] = r;
  }
}

static void
_bspl5d2_Nf(float *f, const float *x, size_t len, const double *parm) {
  float ax, tmp, r;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = AIR_ABS(x[i]);
    BSPL5D2(r, tmp, ax);
    f[i] = r;
  }
}

NrrdKernel
_nrrdKernelBSpline5DD = {
  "bspl5dd",
  BSPL_DECL(5, 2)
};
NrrdKernel *const
nrrdKernelBSpline5DD = &_nrrdKernelBSpline5DD;

/* ------------- order *5* approximate numerical inverse -------------- */

static double
_bspl5_ANI_sup(const double *parm) {
  AIR_UNUSED(parm);
  return 19.5;
}

static double
_bspl5_ANI_int(const double *parm) {
  AIR_UNUSED(parm);
  return 1.0;
}

static double
_bspl5_ANI_kvals[19] = {
  2.842170922021427870236333,
  -1.321729472987239796417307,
  0.5733258709611149890510146,
  -0.2470419274010479815114381,
  0.1063780046404650785440854,
  -0.04580408418467518130037713,
  0.01972212399699206014654736,
  -0.008491860984275658620122180,
  0.003656385950780789716770681,
  -0.001574349495225446217828165,
  0.0006778757185045443332966769,
  -0.0002918757322635763049702028,
  0.0001256725426338698784062181,
  -0.00005410696497728715841372199,
  0.00002328659592249373987497103,
  -0.00001000218170092531503506361,
  4.249940115067599514119408e-6,
  -1.698979738236873388431330e-6,
  4.475539012615912040164139e-7};

#define BSPL5_ANI(ret, tmp, x)                  \
  tmp = AIR_CAST(unsigned int, x+0.5);          \
  if (tmp < 19) {                               \
    ret = _bspl5_ANI_kvals[tmp];                \
  } else {                                      \
    ret = 0.0;                                  \
  }

static double
_bspl5_ANI_1d(double x, const double *parm) {
  double ax, r; int tmp;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL5_ANI(r, tmp, ax);
  return r;
}

static float
_bspl5_ANI_1f(float x, const double *parm) {
  double ax, r; int tmp;
  AIR_UNUSED(parm);

  ax = AIR_ABS(x);
  BSPL5_ANI(r, tmp, ax);
  return AIR_CAST(float, r);
}

static void
_bspl5_ANI_Nd(double *f, const double *x, size_t len, const double *parm) {
  double ax, r; int tmp;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = x[i]; ax = AIR_ABS(ax);
    BSPL5_ANI(r, tmp, ax);
    f[i] = r;
  }
}

static void
_bspl5_ANI_Nf(float *f, const float *x, size_t len, const double *parm) {
  double ax, r; int tmp;
  size_t i;
  AIR_UNUSED(parm);
  
  for (i=0; i<len; i++) {
    ax = x[i]; ax = AIR_ABS(ax);
    BSPL5_ANI(r, tmp, ax);
    f[i] = AIR_CAST(float, r);
  }
}

NrrdKernel
_nrrdKernelBSpline5ApproxInverse = {
  "bspl5ai", 0,
  _bspl5_ANI_sup, _bspl5_ANI_int,
  _bspl5_ANI_1f, _bspl5_ANI_Nf,
  _bspl5_ANI_1d, _bspl5_ANI_Nd
};
NrrdKernel *const
nrrdKernelBSpline5ApproxInverse = &_nrrdKernelBSpline5ApproxInverse;
