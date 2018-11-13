/*
  Teem: Tools to process and visualize scientific data and images              
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

#include "ten.h"
#include "privateTen.h"

#define DOF_NUM 1
#define PARM_NUM 1
static const tenModelParmDesc
parmDesc[] = {
  /* 0 */ {"B0", 0.0, TEN_MODEL_B0_MAX, AIR_FALSE, 0},
};

static void 
simulate(double *dwiSim, const double *parm, const tenExperSpec *espec) {
  unsigned int ii;
  double b0;

  b0 = parm[0];
  for (ii=0; ii<espec->imgNum; ii++) {
    dwiSim[ii] = b0;
  }
  return;
}

static char *
parmSprint(char str[AIR_STRLEN_MED], const double *parm) {
  sprintf(str, "(%g)", parm[0]);
  return str;
}

_TEN_PARM_ALLOC
_TEN_PARM_RAND
_TEN_PARM_STEP
_TEN_PARM_DIST
_TEN_PARM_COPY

static int
parmConvert(double *parmDst, const double *parmSrc,
            const tenModel *modelSrc) {
  int ret;

  ret = 0;
  parmDst[0] = parmSrc[0];
  if (modelSrc == tenModelBall) {
    
  } else if (modelSrc == tenModel1Stick) {

  } else if (modelSrc == tenModelBall1Stick) {

  } else if (modelSrc == tenModel1Cylinder) {

  } else if (modelSrc == tenModel1Tensor2) {

  } else {
    unsigned int ii;
    for (ii=0; ii<PARM_NUM; ii++) {
      parmDst[ii] = AIR_NAN;
    }
    ret = 2;
  }
  return ret;
}

_TEN_SQE
_TEN_SQE_GRAD_CENTDIFF
_TEN_SQE_FIT(tenModelB0)

_TEN_NLL

static void
nllGrad(double *grad, const double *parm,
        const tenExperSpec *espec,
        double *dwiBuff, const double *dwiMeas,
        int rician, double sigma) {

  AIR_UNUSED(grad);
  AIR_UNUSED(parm);
  AIR_UNUSED(espec);
  AIR_UNUSED(dwiBuff);
  AIR_UNUSED(dwiMeas);
  AIR_UNUSED(rician);
  AIR_UNUSED(sigma);
  return;
}

static double
nllFit(double *parm, const tenExperSpec *espec,
       const double *dwiMeas, const double *parmInit,
       int rician, double sigma, int knownB0) {
  unsigned int pp;

  AIR_UNUSED(espec);
  AIR_UNUSED(dwiMeas);
  AIR_UNUSED(rician);
  AIR_UNUSED(sigma);
  AIR_UNUSED(knownB0);
  for (pp=0; pp<PARM_NUM; pp++) {
    parm[pp] = parmInit[pp];
  }
  return 0;
}

tenModel
_tenModelB0 = {
  TEN_MODEL_STR_B0,
  _TEN_MODEL_FIELDS
};
const tenModel *const tenModelB0 = &_tenModelB0;
