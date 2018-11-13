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

#include "unrrdu.h"
#include "privateUnrrdu.h"

/* suffix string that indicates percentile-based min/max */
#define PSUFF "%"

#define INFO "Quantize values to 8, 16, or 32 bits"
char *_unrrdu_quantizeInfoL = 
(INFO ". Input values can be fixed point (e.g. quantizing ushorts down to "
 "uchars) or floating point.  Values are clamped to the min and max before "
 "they are quantized, so there is no risk of getting 255 where you expect 0 "
 "(with unsigned char output, for example).  The min and max can be specified "
 "explicitly (as a regular number), or in terms of percentiles (a number "
 "suffixed with \"" PSUFF "\", no space in between).  This does only linear "
 "quantization.  See also \"unu convert\", \"unu 2op x\", and \"unu 3op "
 "clamp\".\n "
 "* Uses nrrdQuantize");

int
unrrdu_quantizeMain(int argc, const char **argv, char *me, hestParm *hparm) {
  hestOpt *opt = NULL;
  char *out, *err;
  Nrrd *nin, *nout;
  char *_minStr, *_maxStr, *minStr, *maxStr, *mmStr;
  int pret, blind8BitRange, minPerc=AIR_FALSE, maxPerc=AIR_FALSE, *mmPerc;
  unsigned int bits, mmIdx, hbins;
  double min, max, minval, maxval, *mm=NULL;
  NrrdRange *range;
  airArray *mop;

  hestOptAdd(&opt, "b,bits", "bits", airTypeOther, 1, 1, &bits, NULL,
             "Number of bits to quantize down to; determines the type "
             "of the output nrrd:\n "
             "\b\bo \"8\": unsigned char\n "
             "\b\bo \"16\": unsigned short\n "
             "\b\bo \"32\": unsigned int",
             NULL, NULL, &unrrduHestBitsCB);
  hestOptAdd(&opt, "min,minimum", "value", airTypeString, 1, 1,
             &_minStr, "nan",
             "The value to map to zero, given explicitly as a regular number, "
             "*or*, if the number is given with a \"" PSUFF"\" suffix, this "
             "minimum is specified in terms of the percentage of samples in "
             "input that are lower. "
             "\"0" PSUFF "\" means the lowest input value is used, "
             "\"1" PSUFF "\" means that the 1% of the lowest values "
             "are all mapped to zero. "
             "By default (not using this option), the lowest input value is "
             "used.");
  hestOptAdd(&opt, "max,maximum", "value", airTypeString, 1, 1,
             &_maxStr, "nan",
             "The value to map to the highest unsigned integral value, given "
             "explicitly as a regular number, "
             "*or*, if the number is given with "
             "a \"" PSUFF "\" suffix, this maximum is specified in terms of "
             "the percentage of samples in input that are higher. "
             "\"0" PSUFF "\" means the highest input value is used, "
             "which is also the default "
             "behavior (same as not using this option).");
  hestOptAdd(&opt, "hb,bins", "bins", airTypeUInt, 1, 1, &hbins, "5000",
             "number of bins in histogram of values, for determining min "
             "or max by percentiles.  This has to be large enough so that "
             "any errant very high or very low values do not compress the "
             "interesting part of the histogram to an inscrutably small "
             "number of bins.");
  hestOptAdd(&opt, "blind8", "bool", airTypeBool, 1, 1, &blind8BitRange,
             nrrdStateBlind8BitRange ? "true" : "false",
             "if not using \"-min\" or \"-max\", whether to know "
             "the range of 8-bit data blindly (uchar is always [0,255], "
             "signed char is [-128,127])");
  OPT_ADD_NIN(nin, "input nrrd");
  OPT_ADD_NOUT(out, "output nrrd");

  mop = airMopNew();
  airMopAdd(mop, opt, (airMopper)hestOptFree, airMopAlways);

  USAGE(_unrrdu_quantizeInfoL);
  PARSE();
  airMopAdd(mop, opt, (airMopper)hestParseFree, airMopAlways);

  minStr = airStrdup(_minStr);
  airMopAdd(mop, minStr, airFree, airMopAlways);
  maxStr = airStrdup(_maxStr);
  airMopAdd(mop, maxStr, airFree, airMopAlways);

  /* parse min and max */
  min = max = AIR_NAN;
  for (mmIdx=0; mmIdx<=1; mmIdx++) {
    if (0 == mmIdx) {
      mm = &min;
      mmStr = minStr;
      mmPerc = &minPerc;
    } else {
      mm = &max;
      mmStr = maxStr;
      mmPerc = &maxPerc;
    }
    if (airEndsWith(mmStr, PSUFF)) {
      *mmPerc = AIR_TRUE;
      mmStr[strlen(mmStr)-strlen(PSUFF)] = '\0';
    }
    if (1 != airSingleSscanf(mmStr, "%lf", mm)) {
      fprintf(stderr, "%s: couldn't parse \"%s\" as %s\n", me, 
              !mmIdx ? _minStr : _maxStr,
              !mmIdx ? "minimum" : "maximum");
      airMopError(mop);
      return 1;
    }
  }

  if (minPerc || maxPerc) {
    /* have to compute histogram to find real min or max */
    Nrrd *nhist;
    double *hist, sum, total;
    unsigned int hi;

    nhist = nrrdNew();
    airMopAdd(mop, nhist, (airMopper)nrrdNuke, airMopAlways);
    if (nrrdHisto(nhist, nin, NULL, NULL, hbins, nrrdTypeDouble)) {
      airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble making histogram:\n%s", me, err);
      airMopError(mop);
      return 1;
    }
    hist = AIR_CAST(double *, nhist->data);
    total = AIR_CAST(double, nrrdElementNumber(nin));
    if (minPerc) {
      sum = 0;
      for (hi=0; hi<hbins; hi++) {
        sum += hist[hi];
        if (sum >= min*total/100.0) {
          minval = AIR_AFFINE(0, hi, hbins-1,
                              nhist->axis[0].min, nhist->axis[0].max);
          break;
        }
      }
      if (hi == hbins) {
        biffAddf(NRRD, "%s: failed to find lower %g-percentile value",
                 me, min);
        airMopError(mop);
        return 1;
      }
      /* fprintf(stderr, "!%s: %g-%% min = %g\n", me, min, minval); */
    } else {
      minval = min;
    }
    if (maxPerc) {
      sum = 0;
      for (hi=hbins; hi; hi--) {
        sum += hist[hi-1];
        if (sum >= max*total/100.0) {
          maxval = AIR_AFFINE(0, hi-1, hbins-1,
                              nhist->axis[0].min, nhist->axis[0].max);
          break;
        }
      }
      if (!hi) {
        biffAddf(NRRD, "%s: failed to find upper %g-percentile value", me,
                 max);
        return 1;
      }
      /* fprintf(stderr, "!%s: %g-%% max = %g\n", me, max, maxval); */
    } else {
      maxval = max;
    }
  } else {
    minval = min;
    maxval = max;
  }

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);

  /* If the user never specified min or max, they'll be
     AIR_NAN, and nrrdRangeSafeSet will find them, and will do so
     according to blind8BitRange */
  range = nrrdRangeNew(minval, maxval);
  airMopAdd(mop, range, (airMopper)nrrdRangeNix, airMopAlways);
  nrrdRangeSafeSet(range, nin, blind8BitRange);
  if (nrrdQuantize(nout, nin, range, bits)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error quantizing nrrd:\n%s", me, err);
    airMopError(mop);
    return 1;
  }

  SAVE(out, nout, NULL);

  airMopOkay(mop);
  return 0;
}

UNRRDU_CMD(quantize, INFO);
