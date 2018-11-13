/*
 *  _tt_line_backproject_ray_cpu.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, 2019-2013.
 *  CMIC - Centre for Medical Image Computing, UCL, London. (2009-2013)
 *  Aalto University School of Science, Helsinki. (2013) 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_tt_line_backproject_ray_cpu.h"

/*
struct Ray_cpu {
	float_3 o;	// origin
	float_3 d;	// direction
};
int intersectBox_cpu(Ray_cpu r, float_3 boxmin, float_3 boxmax, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float_3 invR = {1.0/r.d.x,1.0/r.d.y,1.0/r.d.z};
    float_3 tbot = {invR.x*(boxmin.x - r.o.x), invR.y*(boxmin.y - r.o.y), invR.z*(boxmin.z - r.o.z)};
    float_3 ttop = {invR.x*(boxmax.x - r.o.x), invR.y*(boxmax.y - r.o.y), invR.z*(boxmax.z - r.o.z)};

    // re-order intersections to find smallest and largest on each axis
    float_3 tmin = fminf(ttop, tbot);
    float_3 tmax = fmaxf(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
    float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

    *tnear = largest_tmin;
    *tfar = smallest_tmax;

    return smallest_tmax > largest_tmin;
}
float_3 mul_cpu(float *M, float_3 v)
{
    float_3 r;
    float_3 t = {M[0], M[1], M[2]};
    r.x = dot(v, t);
    t.x=M[4]; t.y=M[5]; t.z=M[6];
    r.y = dot(v, t);
    t.x=M[8]; t.y=M[9]; t.z=M[10];
    r.z = dot(v, t);
    return r;
}
float_4 mul_cpu(float *M, float_4 v)
{
    float_4 r;
    float_4 t = {M[0], M[1], M[2], M[3]};
    r.x = dot(v, t);
    t.x=M[4]; t.y=M[5]; t.z=M[6]; t.w=M[7];
    r.y = dot(v, t);
    t.x=M[8]; t.y=M[9]; t.z=M[10]; t.w=M[11];
    r.z = dot(v, t);
    r.w = 1.0f;
    return r;}
*/
extern "C" int tt_line_backproject_ray_cpu(float *out_backprojection, float *current_projection, float *invViewMatrix, u_int_2 detectorPixels, float_3 sourcePosition, u_int_3 volumeVoxels, float_3 volumeSize, float t_step, int interpolation)
{return 1;}
/*
{   
    const int    maxSteps = 100000; 
    const float_3 boxMin = {0.0f, 0.0f, 0.0f};
    const float_3 boxMax = {volumeSize.x, volumeSize.y, volumeSize.z};
    float u,v;
    float_4 temp4; 
    float_3 temp3; 
    int hits=0;

    for (int x=1; x<detectorPixels.w; x++)
    {
        for (int y=1; y<detectorPixels.h; y++)
        {
            u = (x / (float) detectorPixels.w);
            v = (y / (float) detectorPixels.h);

            Ray_cpu eyeRay;
            eyeRay.o = sourcePosition;
            //transform and normalize direction vector
            temp4.x=u; temp4.y=v; temp4.z=0.0f; temp4.w=1.0f;
            temp4 = mul_cpu(invViewMatrix, temp4);
            temp3.x = temp4.x; temp3.y = temp4.y; temp3.z = temp4.z; 
            eyeRay.d = normalize(temp3-eyeRay.o); 
            // find intersection with box
            float tnear, tfar;
            int hit = intersectBox_cpu(eyeRay, boxMin, boxMax, &tnear, &tfar);
//            if (!hit) fprintf(stderr,"\n -> %d %d - %f %f %f - %f %f %f - %f %f %f - %f %f %f",x,y,boxMin.x, boxMin.y, boxMin.z, boxMax.x, boxMax.y, boxMax.z, eyeRay.o.x,eyeRay.o.y,eyeRay.o.z, eyeRay.d.x,eyeRay.d.y,eyeRay.d.z );
            if (tnear < 0.0f) tnear = 0.0f;
            if (hit)
            {
                hits++;

                // march along ray from front to back, accumulating
                float  t = tnear;
                float_3 pos = eyeRay.o + eyeRay.d*tnear;
                float_3 step = eyeRay.d*t_step;
                float  bkpr = current_projection[y*detectorPixels.w+x];

                for(int i=0; i<maxSteps; i++)
                {
                    if (interpolation == INTERP_NEAREST)
                    {
                        int px, py, pz;
                        px = (int) pos.x; py = (int) pos.y; pz = (int) pos.z; 
                        int index = pz*volumeVoxels.y*volumeVoxels.x+py*volumeVoxels.x+px;
                        if (px<volumeVoxels.x && px>=0)
                            if (py<volumeVoxels.y && py>=0)
                                if (pz<volumeVoxels.z && pz>=0)
                                    out_backprojection[index] += bkpr;
                    }
                    else if (interpolation == INTERP_TRILINEAR)
                    {
                        int_3 p000, p001, p010, p011, p100, p101, p110, p111;
                        p000.x = (int)floor(pos.x); p000.y = (int)floor(pos.y); p000.z = (int)floor(pos.z); 
                        p001 = p000; p001.x+=1;
                        p010 = p000; p010.y+=1;
                        p011 = p010; p011.x+=1;
                        p100 = p000; p100.z+=1;
                        p101 = p001; p101.z+=1;
                        p110 = p100; p110.y+=1;
                        p111 = p110; p111.x+=1; 
                        float_3 d;
                        d.x = pos.x-p000.x;
                        d.y = pos.y-p000.y;
                        d.z = pos.z-p000.z; 
                        float w000, w001, w010, w011, w100, w101, w110, w111;                      
                        w000 = (1-d.z)*(1-d.y)*(1-d.x)*bkpr;
                        w001 = (1-d.z)*(1-d.y)*( d.x )*bkpr;
                        w010 = (1-d.z)*( d.y )*(1-d.x)*bkpr;
                        w011 = (1-d.z)*( d.y )*( d.x )*bkpr;
                        w100 = ( d.z )*(1-d.y)*(1-d.x)*bkpr;
                        w101 = ( d.z )*(1-d.y)*( d.x )*bkpr;
                        w110 = ( d.z )*( d.y )*(1-d.x)*bkpr;
                        w111 = ( d.z )*( d.y )*( d.x )*bkpr;
                        out_backprojection[p000.z*volumeVoxels.y*volumeVoxels.x+p000.y*volumeVoxels.x+p000.x] += w000;
                        out_backprojection[p001.z*volumeVoxels.y*volumeVoxels.x+p001.y*volumeVoxels.x+p001.x] += w001;
                        out_backprojection[p010.z*volumeVoxels.y*volumeVoxels.x+p010.y*volumeVoxels.x+p010.x] += w010;
                        out_backprojection[p011.z*volumeVoxels.y*volumeVoxels.x+p011.y*volumeVoxels.x+p011.x] += w011;
                        out_backprojection[p100.z*volumeVoxels.y*volumeVoxels.x+p100.y*volumeVoxels.x+p100.x] += w100;
                        out_backprojection[p101.z*volumeVoxels.y*volumeVoxels.x+p101.y*volumeVoxels.x+p101.x] += w101;
                        out_backprojection[p110.z*volumeVoxels.y*volumeVoxels.x+p110.y*volumeVoxels.x+p110.x] += w110;
                        out_backprojection[p111.z*volumeVoxels.y*volumeVoxels.x+p111.y*volumeVoxels.x+p111.x] += w111;
                    }

                    t += t_step;
                    if (t > tfar) break;
                    pos += step;
                }
            }
        }
    }
//    fprintf(stderr,"\nHits: %d",hits);
    return 0;
}

*/

