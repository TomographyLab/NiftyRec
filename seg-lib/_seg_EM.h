#pragma once

#include "_seg_common.h"
#include "_seg_BiasCorrection.h"
#include "_seg_tools.h"
#include "_seg_MRF.h"
#include "_seg_FMM.h"
#include "_seg_Topo.h"
#include <cmath>


class seg_EM
{
protected:

    nifti_image*    inputImage; // pointer to external
    bool    inputImage_status;
    string  FilenameOut;
    int     verbose_level;

    // Size
    int     dimentions;
    int     nx;
    int     ny;
    int     nz;
    int     nt;
    int     nu;
    float     dx;
    float     dy;
    float     dz;
    int     numel;
    int     iter;
    int checkpoint_iter;
    float ratio;
    float reg_factor;
    ImageSize * CurrSizes;


    // SegParameters
    float*  M;
    float*  V;
    float*  Expec;
    float*  ShortPrior;
    int*    Short_2_Long_Indices;
    int*    Long_2_Short_Indices;
    int     numb_classes;
    float   loglik;
    float   oldloglik;
    int     maxIteration;
    bool    aprox;

    // Mask
    nifti_image*    Mask; // pointer to external
    bool    maskImage_status;
    int     numelmasked;

    // Priors Specific
    bool    Priors_status;
    nifti_image*  Priors;

    // MRF Specific
    bool    MRF_status;
    float   MRF_strength;
    float*  MRF;
    float*  MRF_beta;
    float*  MRF_transitionMatrix; // G matrix

    float * Outlierness;
    float * OutliernessUSE;
    bool OutliernessFlag;
    float OutliernessThreshold;
    float Outlierness_ratio;

    // BiasField Specific
    bool    BiasField_status;
    int     BiasField_order;
    float*  BiasField;
    float*  BiasField_coeficients;
    float   BiasField_ratio;
    int     numelbias;

    // LoAd Specific
    int LoAd_phase;
    bool PV_model_status;
    bool SG_deli_status;
    bool Relax_status;
    float  Relax_factor;
    float RelaxGaussKernelSize;
    bool MAP_status;
    float* MAP_M;
    float* MAP_V;

    // Private funcs
    int Create_diagonal_MRF_transitionMatrix();
    int Normalize_T1();
    int Create_CurrSizes();
    int Maximization();
    int Expectation();
    int UpdateMRF();
    int UpdatePriorWeight();
    int UpdateBiasField();
    int UpdateOutlierness();
    int Normalize_Image_and_Priors();
    int Allocate_and_Initialize();
    int Intensity_Based_Inisitalization_of_Means();

public:
    seg_EM(int numb_classes,int NumbMultiSpec,int NumbTimePoints);
    ~seg_EM();
    int SetInputImage(nifti_image *);
    int SetMaskImage(nifti_image *);
    int SetPriorImage(nifti_image *);
    int SetFilenameOut(const char *);
    int SetMAP(float *M, float* V);
    int SetRegValue(float reg);
    int Turn_Relaxation_ON(float relax_factor,float relax_gauss_kernel);
    int Turn_MRF_ON(float MRF_strenght);
    int OutliernessON(float OutliernessThreshold, float ratio);
    int Turn_BiasField_ON(int BiasField_order,float ratiothresh);
    int SetLoAd(float RelaxFactor,bool Relax_status,bool PV_model_status,bool SG_deli_status);
    int SetMaximalIterationNumber(unsigned int numberiter);
    int SetAprox(bool aproxval);
    int SetVerbose(unsigned int verblevel);
    int FreeInputImage(void);
    int FreePriorImage(void);
    int FreeMaskImage(void);
    int FreeBiasImage(void);

    int CheckParameters_EM();
    int Initisalise_EM();
    int Run_EM();
    int Init_EM();
    int Step_EM(int n);
    int Step_EM_Gaussian();
    int Step_EM_Expectation();
    int Step_EM_Maximization();
    int Step_EM_MRF();
    int Step_EM_BiasField();
    int Step_EM_PriorWeight();
    float * GetExpec();
    int ChangeInputImage(SegPrecisionTYPE* image_data);
    nifti_image * GetInputImage();
    int get_size_x();
    int get_size_y();
    int get_size_z();
    int get_n_images();
    int get_n_classes();
    int SetExpec(SegPrecisionTYPE * expec_data);
    float GetLoglik();
    float * GetMean();
    float * GetVariance();
    float GetRegValue();
    int SetMean(float *mean);
    int SetVariance(float *variance);

    nifti_image *GetMask();
    nifti_image *GetPriorImage();
    nifti_image *GetResult();
    nifti_image *GetResultNeonate();
    nifti_image *GetBiasCorrected(const char * filename);
    nifti_image *GetOutlierness(const char * filename);
};

