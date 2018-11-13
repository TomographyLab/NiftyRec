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

#include "gage.h"
#include "privateGage.h"

/* these functions don't necessarily belong in gage, but we're putting
   them here for the time being.  Being in gage means that vprobe and
   pprobe don't need extra libraries to find them. */

#define BT 2.526917043979558
#define AA 0.629078014852877

#define TAU_OF_TEE(tee) (tee < BT ? AA*sqrt(tee) : 0.5365 + log(tee)/2)
/* is it surprising that the transition is at tau = 1 ? */
#define TEE_OF_TAU(tau) (tau < 1 ? tau*tau/(AA*AA) : exp(2*(tau - 0.5365)))

double
gageTauOfTee(double tee) {
  double tau;

  tau = TAU_OF_TEE(tee);
  return tau;
}

double
gageTeeOfTau(double tau) {
  double tee;

  tee = TEE_OF_TAU(tau);
  return tee;
}

double
gageSigOfTau(double tau) {
  
  return sqrt(TEE_OF_TAU(tau));
}

double
gageTauOfSig(double sig) {

  return TAU_OF_TEE(sig*sig);
}

#undef BT
#undef AA
#undef TAU_OF_TEE
#undef TEE_OF_TAU

double
gageStackWtoI(gageContext *ctx, double swrl, int *outside) {
  double si;

  if (ctx && ctx->parm.stackUse && outside) {
    unsigned int sidx;
    if (swrl < ctx->stackPos[0]) {
      /* we'll extrapolate from stackPos[0] and [1] */
      sidx = 0;
      *outside = AIR_TRUE;
    } else if (swrl > ctx->stackPos[ctx->pvlNum-2]) {
      /* extrapolate from stackPos[ctx->pvlNum-3] and [ctx->pvlNum-2];
         gageStackPerVolumeAttach ensures that we there are at least two
         blurrings pvls & one base pvl ==> pvlNum >= 3 ==> pvlNum-3 >= 0 */
      sidx = ctx->pvlNum-3;
      *outside = AIR_TRUE;
    } else {
      /* HEY: stupid linear search */
      for (sidx=0; sidx<ctx->pvlNum-2; sidx++) {
        if (AIR_IN_CL(ctx->stackPos[sidx], swrl, ctx->stackPos[sidx+1])) {
          break;
        }
      }
      if (sidx == ctx->pvlNum-2) {
        /* search failure */
        *outside = AIR_FALSE;
        return AIR_NAN;
      }
      *outside = AIR_FALSE;
    }
    si = AIR_AFFINE(ctx->stackPos[sidx], swrl, ctx->stackPos[sidx+1],
                    sidx, sidx+1);
  } else {
    si = AIR_NAN;
  }
  return si;
}

double
gageStackItoW(gageContext *ctx, double si, int *outside) {
  unsigned int sidx;
  double swrl, sfrac;

  if (ctx && ctx->parm.stackUse && outside) {
    if (si < 0) {
      sidx = 0;
      *outside = AIR_TRUE;
    } else if (si > ctx->pvlNum-2) {
      sidx = ctx->pvlNum-3;
      *outside = AIR_TRUE;
    } else {
      sidx = AIR_CAST(unsigned int, si);
      *outside = AIR_FALSE;
    }
    sfrac = si - sidx;
    swrl = AIR_AFFINE(0, sfrac, 1, ctx->stackPos[sidx], ctx->stackPos[sidx+1]);
    /*
    fprintf(stderr, "!%s: si %g (%u) --> %u + %g --> [%g,%g] -> %g\n", me,
            si, ctx->pvlNum, sidx, sfrac, 
            ctx->stackPos[sidx], ctx->stackPos[sidx+1], swrl);
    */
  } else {
    swrl = AIR_NAN;
  }
  return swrl;
}

/*
** this is a little messy: the given pvlStack array has to be allocated
** by the caller to hold blNum gagePerVolume pointers, BUT, the values
** of pvlStack[i] shouldn't be set to anything: as with gagePerVolumeNew(),
** gage allocates the pervolume itself.
*/
int
gageStackPerVolumeNew(gageContext *ctx,
                      gagePerVolume **pvlStack,
                      const Nrrd *const *nblur, unsigned int blNum,
                      const gageKind *kind) {
  static const char me[]="gageStackPerVolumeNew";
  unsigned int blIdx;

  if (!( ctx && pvlStack && nblur && kind )) {
    biffAddf(GAGE, "%s: got NULL pointer", me);
    return 1;
  }
  if (!blNum) {
    biffAddf(GAGE, "%s: need non-zero num", me);
    return 1;
  }

  for (blIdx=0; blIdx<blNum; blIdx++) {
    if (!( pvlStack[blIdx] = gagePerVolumeNew(ctx, nblur[blIdx], kind) )) {
      biffAddf(GAGE, "%s: on pvl %u of %u", me, blIdx, blNum);
      return 1;
    }
  }

  return 0;
}

/*
** the "base" pvl is the LAST pvl, ctx->pvl[pvlNum-1]
*/
int
gageStackPerVolumeAttach(gageContext *ctx, gagePerVolume *pvlBase,
                         gagePerVolume **pvlStack, const double *stackPos,
                         unsigned int blNum) {
  static const char me[]="gageStackPerVolumeAttach";
  unsigned int blIdx;

  if (!(ctx && pvlBase && pvlStack && stackPos)) { 
    biffAddf(GAGE, "%s: got NULL pointer %p %p %p %p", me,
             ctx, pvlBase, pvlStack, stackPos);
    return 1;
  }
  if (!( blNum >= 2 )) {
    /* this constraint is important for the logic of stack reconstruction:
       minimum number of node-centered samples is 2, and the number of
       pvls has to be at least 3 (two blurrings + one base pvl) */
    biffAddf(GAGE, "%s: need at least two samples along stack", me);
    return 1;
  }
  if (ctx->pvlNum) {
    biffAddf(GAGE, "%s: can't have pre-existing volumes (%u) "
             "prior to stack attachment", me, ctx->pvlNum);
    return 1;
  }
  for (blIdx=0; blIdx<blNum; blIdx++) {
    if (!AIR_EXISTS(stackPos[blIdx])) {
      biffAddf(GAGE, "%s: stackPos[%u] = %g doesn't exist", me, blIdx, 
               stackPos[blIdx]);
      return 1;
    }
    if (blIdx < blNum-1) {
      if (!( stackPos[blIdx] < stackPos[blIdx+1] )) {
        biffAddf(GAGE, "%s: stackPos[%u] = %g not < stackPos[%u] = %g", me,
                 blIdx, stackPos[blIdx], blIdx+1, stackPos[blIdx+1]);
        return 1;
      }
    }
  }

  /* the base volume is LAST, after all the stack samples */
  for (blIdx=0; blIdx<blNum; blIdx++) {
    if (gagePerVolumeAttach(ctx, pvlStack[blIdx])) {
      biffAddf(GAGE, "%s: on pvl %u of %u", me, blIdx, blNum);
      return 1;
    }
  }
  if (gagePerVolumeAttach(ctx, pvlBase)) {
    biffAddf(GAGE, "%s: on base pvl", me);
    return 1;
  }
  
  airFree(ctx->stackPos);
  airFree(ctx->stackFsl);
  airFree(ctx->stackFw);
  ctx->stackPos = AIR_CALLOC(blNum, double);
  ctx->stackFsl = AIR_CALLOC(blNum, double);
  ctx->stackFw = AIR_CALLOC(blNum, double);
  if (!( ctx->stackPos && ctx->stackFsl && ctx->stackFw )) {
    biffAddf(GAGE, "%s: couldn't allocate stack buffers (%p %p %p)", me,
             AIR_CAST(void *, ctx->stackPos),
             AIR_CAST(void *, ctx->stackFsl),
             AIR_CAST(void *, ctx->stackFw));
    return 1;
  }
  for (blIdx=0; blIdx<blNum; blIdx++) {
    ctx->stackPos[blIdx] = stackPos[blIdx];
  }

  return 0;
}

/*
** _gageStackBaseIv3Fill
**
** after the individual iv3's in the stack have been filled, this does
** the across-stack filtering to fill the iv3 of pvl[pvlNum-1] (the
** "base" pvl) 
*/
int
_gageStackBaseIv3Fill(gageContext *ctx) {
  static const char me[]="_gageStackBaseIv3Fill";
  unsigned int fd, pvlIdx, cacheIdx, cacheLen, baseIdx, valLen;

  fd = 2*ctx->radius;
  /* the "base" pvl is the LAST pvl */
  baseIdx = ctx->pvlNum - 1; 
  cacheLen = fd*fd*fd*ctx->pvl[0]->kind->valLen;
  if (ctx->verbose > 2) {
    fprintf(stderr, "%s: cacheLen = %u\n", me, cacheLen);
  }
  if (nrrdKernelHermiteFlag == ctx->ksp[gageKernelStack]->kernel) {
    unsigned int xi, yi, zi, blurIdx, valIdx, fdd;
    double xx, *iv30, *iv31, sigma0, sigma1;
    
    fdd = fd*fd;
    /* initialize the output iv3 to all zeros, since we won't be
       usefully setting the values on the boundary (the boundary which
       is required in the rest of the stack's iv3s in order to do the
       laplacian-based spline recon), and we can't have any
       non-existent values creeping in.  We shouldn't need to do any
       kind of nrrdBoundaryBleed thing here, because the kernel
       weights really should be zero on the boundary. */
    for (cacheIdx=0; cacheIdx<cacheLen; cacheIdx++) {
      ctx->pvl[baseIdx]->iv3[cacheIdx] = 0;
    }

    /* find the interval in the pre-blurred volumes containing the
       desired scale location */
    for (pvlIdx=0; pvlIdx<ctx->pvlNum-1; pvlIdx++) {
      if (ctx->stackFw[pvlIdx]) {
        /* has to be non-zero somewhere, since _gageLocationSet()
           gives an error if there aren't non-zero stackFw[i] */
        break;
      }
    }
    /* so no way that pvlIdx == pvlNum-1 */
    if (pvlIdx == ctx->pvlNum-2) {
      /* pvlNum-2 is pvl index of last pre-blurred volume */
      /* gageStackPerVolumeAttach() enforces getting at least two 
         pre-blurred volumes --> pvlNum >= 3 --> blurIdx >= 0 */
      blurIdx = pvlIdx-1;
      xx = 1;
    } else {
      blurIdx = pvlIdx;
      /* by design, the hermite non-kernel generates the same values as
         the tent kernel (with scale forced == 1), so we can use those
         to control the interpolation */
      xx = 1 - ctx->stackFw[pvlIdx];
    }
    iv30 = ctx->pvl[blurIdx]->iv3;
    iv31 = ctx->pvl[blurIdx+1]->iv3;
    sigma0 = ctx->stackPos[blurIdx];
    sigma1 = ctx->stackPos[blurIdx+1];
    valLen = ctx->pvl[baseIdx]->kind->valLen;
    for (valIdx=0; valIdx<valLen; valIdx++) {
      unsigned iii;
      double val0, val1, drv0, drv1, lapl0, lapl1, aa, bb, cc, dd;
      for (zi=1; zi<fd-1; zi++) {
        for (yi=1; yi<fd-1; yi++) {
          for (xi=1; xi<fd-1; xi++) {
            /* note that iv3 axis ordering is x, y, z, tuple */
            iii = xi + fd*(yi + fd*(zi + fd*valIdx));
            val0 = iv30[iii];
            val1 = iv31[iii];
            lapl0 = (iv30[iii + 1]   + iv30[iii - 1] +
                     iv30[iii + fd]  + iv30[iii - fd] +
                     iv30[iii + fdd] + iv30[iii - fdd] - 6*val0);
            lapl1 = (iv31[iii + 1]   + iv31[iii - 1] +
                     iv31[iii + fd]  + iv31[iii - fd] +
                     iv31[iii + fdd] + iv31[iii - fdd] - 6*val1);
            /* the (sigma1 - sigma0) factor is needed to convert the
               derivative with respect to sigma (sigma*lapl) into the
               derivative with respect to xx */
            drv0 = sigma0*lapl0*(sigma1 - sigma0);
            drv1 = sigma1*lapl1*(sigma1 - sigma0);
            /* Hermite spline coefficients, thanks Mathematica */
            aa = drv0 + drv1 + 2*val0 - 2*val1;
            bb = -2*drv0 - drv1 - 3*val0 + 3*val1;
            cc = drv0;
            dd = val0;
            ctx->pvl[baseIdx]->iv3[iii] = dd + xx*(cc + xx*(bb + aa*xx));
          }
        }
      }
    }
  } else {
    /* we're doing simple convolution-based recon on the stack */
    /* NOTE we are treating the 4D fd*fd*fd*valLen iv3 as a big 1-D array */
    double wght, val;
    for (cacheIdx=0; cacheIdx<cacheLen; cacheIdx++) {
      val = 0;
      for (pvlIdx=0; pvlIdx<ctx->pvlNum-1; pvlIdx++) {
        wght = ctx->stackFw[pvlIdx];
        val += (wght
                ? wght*ctx->pvl[pvlIdx]->iv3[cacheIdx]
                : 0);
      }
      ctx->pvl[baseIdx]->iv3[cacheIdx] = val;
    }
  }
  return 0;
}

/*
******** gageStackProbe()
*/
int
gageStackProbe(gageContext *ctx,
               double xi, double yi, double zi, double stackIdx) {
  static const char me[]="gageStackProbe";

  if (!ctx) {
    return 1;
  }
  if (!ctx->parm.stackUse) {
    sprintf(ctx->errStr, "%s: can't probe stack without parm.stackUse", me);
    ctx->errNum = 1;
    return 1;
  }
  return _gageProbe(ctx, xi, yi, zi, stackIdx);
}

int
gageStackProbeSpace(gageContext *ctx,
                    double xx, double yy, double zz, double ss,
                    int indexSpace, int clamp) {
  static const char me[]="gageStackProbeSpace";

  if (!ctx) {
    return 1;
  }
  if (!ctx->parm.stackUse) {
    sprintf(ctx->errStr, "%s: can't probe stack without parm.stackUse", me);
    ctx->errNum = 1;
    return 1;
  }
  return _gageProbeSpace(ctx, xx, yy, zz, ss, indexSpace, clamp);
}

