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

#define INFO "Simulate DW images from a field of models"
char *_tend_msimInfoL =
  (INFO
   ".  The output will be in the same form as the input to \"tend estim\". "
   "The B-matrices (\"-B\") can be the output from \"tend bmat\", or the "
   "gradients can be given directly (\"-g\"); one of these is required. "
   "Note that the input tensor field (\"-i\") is the basis of the output "
   "per-axis fields and image orientation.  NOTE: this includes the "
   "measurement frame used in the input tensor field, which implies that "
   "the given gradients or B-matrices are already expressed in that "
   "measurement frame. ");

int
tend_msimMain(int argc, const char **argv, char *me, hestParm *hparm) {
  int pret;
  hestOpt *hopt = NULL;
  char *perr, *err;
  airArray *mop;

  tenExperSpec *espec;
  const tenModel *model;
  int E, seed, keyValueSet, outType, plusB0;
  Nrrd *nin, *nT2, *_ngrad, *ngrad, *nout;
  char *outS;
  double bval, sigma;

  /* maybe this can go in tend.c, but for some reason its explicitly
     set to AIR_FALSE there */
  hparm->elideSingleOtherDefault = AIR_TRUE;

  hestOptAdd(&hopt, "sigma", "sigma", airTypeDouble, 1, 1, &sigma, "0.0",
             "Rician noise parameter");
  hestOptAdd(&hopt, "seed", "seed", airTypeInt, 1, 1, &seed, "42",
             "seed value for RNG which creates noise");
  hestOptAdd(&hopt, "g", "grad list", airTypeOther, 1, 1, &_ngrad, NULL,
             "gradient list, one row per diffusion-weighted image", 
             NULL, NULL, nrrdHestNrrd);
  hestOptAdd(&hopt, "r", "reference field", airTypeOther, 1, 1, &nT2, NULL,
             "reference anatomical scan, with no diffusion weighting",
             NULL, NULL, nrrdHestNrrd);
  hestOptAdd(&hopt, "i", "model field", airTypeOther, 1, 1, &nin, "-",
             "input model field", NULL, NULL, nrrdHestNrrd);
  hestOptAdd(&hopt, "b", "b", airTypeDouble, 1, 1, &bval, "1000",
             "b value for simulated scan");
  hestOptAdd(&hopt, "kvp", NULL, airTypeInt, 0, 0, &keyValueSet, NULL,
             "generate key/value pairs in the NRRD header corresponding "
             "to the input b-value and gradients.");
  hestOptAdd(&hopt, "t", "type", airTypeEnum, 1, 1, &outType, "float",
             "output type of DWIs", NULL, nrrdType);
  hestOptAdd(&hopt, "o", "nout", airTypeString, 1, 1, &outS, "-",
             "output dwis");

  mop = airMopNew();
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  USAGE(_tend_msimInfoL);
  PARSE();
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);
  espec = tenExperSpecNew();
  airMopAdd(mop, espec, (airMopper)tenExperSpecNix, airMopAlways);

  airSrandMT(seed);
  if (nrrdTypeDouble == _ngrad->type) {
    ngrad = _ngrad;
  } else {
    ngrad = nrrdNew();
    airMopAdd(mop, ngrad, (airMopper)nrrdNuke, airMopAlways);
    if (nrrdConvert(ngrad, _ngrad, nrrdTypeDouble)) {
      airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble converting grads to %s:\n%s\n", me,
              airEnumStr(nrrdType, nrrdTypeDouble), err);
      airMopError(mop); return 1;
    }
  }
  E = 0;
  if (!E) E |= tenGradientCheck(ngrad, nrrdTypeDouble, 1);
  if (!E) E |= tenExperSpecGradSingleBValSet(espec, ngrad->axis[1].size,
                                             bval, ngrad->data);
  if (!E) E |= tenModelFromAxisLearn(&model, &plusB0, nin->axis + 0);
  /* why not do something with plusB0? */
  if (!E) E |= tenModelSimulate(nout, outType, espec,
                                model, nT2, nin, keyValueSet);
  if (E) {
    airMopAdd(mop, err=biffGetDone(TEN), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble:\n%s\n", me, err);
    airMopError(mop); return 1;
  }
  if (nrrdSave(outS, nout, NULL)) {
    airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble writing:\n%s\n", me, err);
    airMopError(mop); return 1;
  }

  airMopOkay(mop);
  return 0;
}
/* TEND_CMD(msim, INFO); */
unrrduCmd tend_msimCmd = { "msim", INFO, tend_msimMain };
