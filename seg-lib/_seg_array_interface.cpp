
#include "_seg_array_interface.h"
#include "_seg_EM.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <string.h>

using namespace std;



typedef struct
{
    bool initialised;
    SEG_PARAM *segment_param;
    seg_EM *SEG;
}
Segmentation;

static Segmentation segmentation;
Segmentation& GlobSegmentation()
{
    return segmentation;
}


bool NiftySeg_initialised()
{
    if (!GlobSegmentation().initialised)
    {
        fprintf(stderr,"NiftySeg not initialised.\n");
    }
    return GlobSegmentation().initialised;
}

/*
void Merge_Priors(nifti_image * Priors, nifti_image ** Priors_temp, SEG_PARAM * segment_param)
{
    long img_size= Priors->nx * Priors->ny * Priors->nz;
    PrecisionTYPE * Prior_ptr_start = static_cast<PrecisionTYPE *>(Priors->data);
    for(long cl=0;cl<(segment_param->numb_classes);cl++){
//printf("Merging prior %d/%d ..\n",cl,(int)segment_param->numb_classes);
        PrecisionTYPE * Prior_tmp_ptr = static_cast<PrecisionTYPE *>(Priors_temp[cl]->data);
        PrecisionTYPE * Prior_ptr = &Prior_ptr_start[cl*img_size];
        //memcpy((void*)Prior_ptr , (const void*)Prior_tmp_ptr , img_size*sizeof(PrecisionTYPE));
        for(long i=0; i<img_size;i++){
            (*Prior_ptr)=(*Prior_tmp_ptr);
//            (*Prior_ptr)=0;
            Prior_ptr++;
            Prior_tmp_ptr++;
        }
    }
    return;
} */

void Merge_Priors(nifti_image * Priors, nifti_image ** Priors_temp, SEG_PARAM * segment_param)
{
    long img_size= Priors->nx * Priors->ny * Priors->nz;
    PrecisionTYPE * Prior_ptr_start = static_cast<PrecisionTYPE *>(Priors->data);
    for(long cl=0;cl<(segment_param->numb_classes);cl++){
        PrecisionTYPE * Prior_tmp_ptr = static_cast<PrecisionTYPE *>(Priors_temp[cl]->data);
        PrecisionTYPE * Prior_ptr = &Prior_ptr_start[cl*img_size];
        for(long i=0; i<img_size;i++){
            (*Prior_ptr)=(*Prior_tmp_ptr);
            Prior_ptr++;
            Prior_tmp_ptr++;
        }
    }
    return;
}




extern "C" int seg_array_initialise_fromfile(int numb_classes, char* filename, char *mask_filename)
{
    /* default parameters */
    GlobSegmentation().segment_param = new SEG_PARAM [1]();
    GlobSegmentation().segment_param->maxIteration=100;
    GlobSegmentation().segment_param->flag_T1=0;
    GlobSegmentation().segment_param->flag_out=0;
    GlobSegmentation().segment_param->flag_mask=1;
    GlobSegmentation().segment_param->flag_MRF=1;
    GlobSegmentation().segment_param->flag_Bias=1;
    GlobSegmentation().segment_param->flag_SG_deli=1;
    GlobSegmentation().segment_param->flag_bc_out=0;
    GlobSegmentation().segment_param->relax_factor=0;
    GlobSegmentation().segment_param->relax_gauss_kernel=0;
    GlobSegmentation().segment_param->flag_PV_model=1;
    GlobSegmentation().segment_param->verbose_level=0;
    GlobSegmentation().segment_param->flag_manual_priors=0;
    GlobSegmentation().segment_param->bias_order=3;
    GlobSegmentation().segment_param->MRF_strength=0.25f;
    GlobSegmentation().segment_param->Bias_threshold=0;
    GlobSegmentation().segment_param->numb_classes=0;
    GlobSegmentation().segment_param->aprox=false;
    GlobSegmentation().segment_param->flag_Outlierness=0;
    GlobSegmentation().segment_param->OutliernessThreshold=0;
    GlobSegmentation().segment_param->flag_out_outlier=0;
    float OutliernessRatio=0.01;

    /* set parameters given as arguments */
    GlobSegmentation().segment_param->numb_classes=numb_classes;
    GlobSegmentation().segment_param->filename_T1 = filename;
    GlobSegmentation().segment_param->flag_T1=1;
    GlobSegmentation().segment_param->filename_out = "seg_result.nii";
    GlobSegmentation().segment_param->flag_out=1;

    /* READ image */
//    nifti_image *T1;
//    int dim[8];
//    dim[0] = 3;
//    dim[1] = size_x;
//    dim[2] = size_y;
//    dim[3] = size_z;
//    dim[4] = 1;
//    dim[5] = 1;
//    dim[6] = 1;
//    dim[7] = 1;
//    T1 = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);  // not templated
//    seg_changeDatatype<PrecisionTYPE>(T1);
//    T1->data = (void*) calloc(size_x*size_y*size_z*n_dimensions,sizeof(PrecisionTYPE));
//    memcpy((void*) T1->data, (void*) image_data, size_x*size_y*size_z*n_dimensions*sizeof(PrecisionTYPE));

    nifti_image * T1=nifti_image_read(GlobSegmentation().segment_param->filename_T1,true);
    if(T1 == NULL){
        fprintf(stderr,"* Error when reading the T1 image: %s\n",GlobSegmentation().segment_param->filename_T1);
        return 1;
    }
    seg_changeDatatype<PrecisionTYPE>(T1);


    /* READ MASK image */
    GlobSegmentation().segment_param->filename_mask = mask_filename;
    GlobSegmentation().segment_param->flag_mask=1;

    nifti_image * Mask=NULL;
    Mask = nifti_image_read(GlobSegmentation().segment_param->filename_mask,true);
    seg_changeDatatype<PrecisionTYPE>(Mask);
    if(Mask == NULL)
    {
        fprintf(stderr,"* Error when reading the mask image: %s\n",GlobSegmentation().segment_param->filename_mask);
        return 1;
    }

    /* Initialised EM classifier */ 
    static seg_EM SEG(GlobSegmentation().segment_param->numb_classes,T1->dim[4],1);
    GlobSegmentation().SEG = &SEG;
    GlobSegmentation().SEG->SetInputImage(T1);

    if(GlobSegmentation().segment_param->flag_mask) {
        GlobSegmentation().SEG->SetMaskImage(Mask);
    }


    if(GlobSegmentation().segment_param->flag_Outlierness){
        GlobSegmentation().SEG->OutliernessON(GlobSegmentation().segment_param->OutliernessThreshold,OutliernessRatio);
    } 

    GlobSegmentation().SEG->SetVerbose(GlobSegmentation().segment_param->verbose_level); 
    GlobSegmentation().SEG->SetFilenameOut(GlobSegmentation().segment_param->filename_out); 

    if(GlobSegmentation().segment_param->flag_Bias) 
        GlobSegmentation().SEG->Turn_BiasField_ON(GlobSegmentation().segment_param->bias_order,GlobSegmentation().segment_param->Bias_threshold); 
    if(GlobSegmentation().segment_param->flag_MRF) 
        GlobSegmentation().SEG->Turn_MRF_ON(GlobSegmentation().segment_param->MRF_strength); 
    if(GlobSegmentation().segment_param->relax_factor>0) 
        GlobSegmentation().SEG->Turn_Relaxation_ON(GlobSegmentation().segment_param->relax_factor,GlobSegmentation().segment_param->relax_gauss_kernel); 
    if(GlobSegmentation().segment_param->flag_MAP) 
        GlobSegmentation().SEG->SetMAP(GlobSegmentation().segment_param->MAP_M,GlobSegmentation().segment_param->MAP_V); 
    GlobSegmentation().SEG->SetAprox(GlobSegmentation().segment_param->aprox); 
    GlobSegmentation().SEG->SetMaximalIterationNumber(GlobSegmentation().segment_param->maxIteration); 

    GlobSegmentation().segment_param->filename_priors= (char **) calloc(GlobSegmentation().segment_param->numb_classes,sizeof(char *)); 

    GlobSegmentation().SEG->Init_EM(); 

    GlobSegmentation().initialised=1; 
    fprintf(stderr,"NiftySeg Initialised.\n"); 
    return 0; 
}


extern "C" int seg_array_set_MRF_strength(double MRF_strength)
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().segment_param->flag_MRF = 1;
    GlobSegmentation().segment_param->MRF_strength = MRF_strength;
    if(GlobSegmentation().segment_param->MRF_strength<=0)
        GlobSegmentation().segment_param->flag_MRF=0;
    GlobSegmentation().SEG->Turn_MRF_ON(GlobSegmentation().segment_param->MRF_strength);
    return 0;
}

extern "C" int seg_array_get_MRF_strength(double *MRF_strength)
{
    if (!NiftySeg_initialised())
        return 1;

    MRF_strength[0] = GlobSegmentation().segment_param->MRF_strength;
    return 0;
}

extern "C" int seg_array_set_regularisation_covariance(double parameter)
{
    if (!NiftySeg_initialised())
        return 1;
    GlobSegmentation().SEG->SetRegValue(parameter);
    return 0;
}

extern "C" int seg_array_get_regularisation_covariance(double *parameter)
{
    if (!NiftySeg_initialised())
        return 1;

    parameter[0] = GlobSegmentation().SEG->GetRegValue();
    return 0;
}

extern "C" int seg_array_set_biasfield_parameters(unsigned int order, double threshold)
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().segment_param->bias_order = order;
    GlobSegmentation().segment_param->Bias_threshold = threshold;
    GlobSegmentation().SEG->Turn_BiasField_ON(GlobSegmentation().segment_param->bias_order,GlobSegmentation().segment_param->Bias_threshold);
    return 0;
}

extern "C" int seg_array_get_biasfield_parameters(int *enabled, unsigned int *order, double *threshold)
{
    if (!NiftySeg_initialised())
        return 1;

    *order = GlobSegmentation().segment_param->bias_order;
    *threshold = GlobSegmentation().segment_param->Bias_threshold;
    if (GlobSegmentation().segment_param->bias_order == 0)
        *enabled = 0;
    else
        *enabled = 1;
    return 0;
}

extern "C" int seg_array_get_biasfield_size(int* size_x, int* size_y, int* size_z)
{
    if (!NiftySeg_initialised())
        return 1;

    *size_x = GlobSegmentation().SEG->get_size_x();
    *size_y = GlobSegmentation().SEG->get_size_y();
    *size_z = GlobSegmentation().SEG->get_size_z();
    return 0;
}

extern "C" int seg_array_get_biasfield(PrecisionTYPE* biasfield_data)
{
    if (!NiftySeg_initialised())
        return 1;

    nifti_image * BiasField = GlobSegmentation().SEG->GetBiasCorrected("bias_corrected.nii");

    int size_x = GlobSegmentation().SEG->get_size_x();
    int size_y = GlobSegmentation().SEG->get_size_y();
    int size_z = GlobSegmentation().SEG->get_size_z();

    memcpy((void*)biasfield_data, (void*)BiasField->data, size_x*size_y*size_z*sizeof(PrecisionTYPE));
    
    return 0;
}

extern "C" int seg_array_save_biasfield(char *filename)
{
    if (!NiftySeg_initialised())
        return 1;

    nifti_image * BiasFieldCorrected=NULL;
    BiasFieldCorrected = GlobSegmentation().SEG->GetBiasCorrected(filename);
    nifti_image_write(BiasFieldCorrected);

    return 0;
}



extern "C" int seg_array_set_image(int size_x, int size_y, int size_z, int n_images, PrecisionTYPE* image_data)
{
    if (!NiftySeg_initialised())
        return 1;

    return 0;
}


extern "C" int seg_array_step(int n)
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().SEG->Step_EM(n);
    return 0;
}

extern "C" int seg_array_step_Gaussian()
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().SEG->Step_EM_Gaussian();
    return 0;
}

extern "C" int seg_array_step_Expectation()
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().SEG->Step_EM_Expectation();
    return 0;
}

extern "C" int seg_array_step_Maximization()
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().SEG->Step_EM_Maximization();
    return 0;
}

extern "C" int seg_array_step_MRF()
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().SEG->Step_EM_MRF();
    return 0;
}
extern "C" int seg_array_step_BiasField()
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().SEG->Step_EM_BiasField();
    return 0;
}
extern "C" int seg_array_step_PriorWeight()
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().SEG->Step_EM_PriorWeight();
    return 0;
}




extern "C" int seg_array_get_image(PrecisionTYPE* image_data)
{
    if (!NiftySeg_initialised())
        return 1;

    nifti_image * Result = GlobSegmentation().SEG->GetInputImage();

    int size_x = GlobSegmentation().SEG->get_size_x();
    int size_y = GlobSegmentation().SEG->get_size_y();
    int size_z = GlobSegmentation().SEG->get_size_z();
    int n_images = GlobSegmentation().SEG->get_n_images();

    memcpy((void*)image_data, (void*)Result->data, size_x*size_y*size_z*n_images*sizeof(PrecisionTYPE));

    return 0;
}

extern "C" int seg_array_get_image_size(int* size_x, int* size_y, int* size_z, int* n_images)
{
    if (!NiftySeg_initialised())
        return 1;

    *size_x = GlobSegmentation().SEG->get_size_x();
    *size_y = GlobSegmentation().SEG->get_size_y();
    *size_z = GlobSegmentation().SEG->get_size_z();
    *n_images = GlobSegmentation().SEG->get_n_images();
    return 0;
}


extern "C" int seg_array_set_input_image(PrecisionTYPE* image_data)
{
    int status=0;
    if (!NiftySeg_initialised())
        return 1;

    status = GlobSegmentation().SEG->ChangeInputImage(image_data);
    return status;
}


extern "C" int seg_array_get_segmentation(PrecisionTYPE* segmentation_data)
{
    if (!NiftySeg_initialised())
        return 1;

    nifti_image * Result = GlobSegmentation().SEG->GetResult();

    int size_x = GlobSegmentation().SEG->get_size_x();
    int size_y = GlobSegmentation().SEG->get_size_y();
    int size_z = GlobSegmentation().SEG->get_size_z();
    int n_classes = GlobSegmentation().SEG->get_n_classes();

    memcpy((void*)segmentation_data, (void*)Result->data, size_x*size_y*size_z*n_classes*sizeof(PrecisionTYPE));

    nifti_image_free(Result);       
    return 0;
}

extern "C" int seg_array_get_segmentation_size(int* size_x, int* size_y, int* size_z, int* n_classes)
{
    if (!NiftySeg_initialised())
        return 1;

    *size_x = GlobSegmentation().SEG->get_size_x();
    *size_y = GlobSegmentation().SEG->get_size_y();
    *size_z = GlobSegmentation().SEG->get_size_z();
    *n_classes = GlobSegmentation().SEG->get_n_classes();
    return 0;
}

extern "C" int seg_array_save_segmentation(char *filename)
{ 
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().SEG->SetFilenameOut(filename);
    nifti_image * Result = GlobSegmentation().SEG->GetResult();
//    Result->fname = filename;
    nifti_image_write(Result);
    nifti_image_free(Result);        
    cout << "Segmentation saved in: " << filename << endl;
    return 0;
}


extern "C" int seg_array_set_segmentation(PrecisionTYPE* segmentation_data)
{
    int status=0;
    if (!NiftySeg_initialised())
        return 1;
    status = GlobSegmentation().SEG->SetExpec(segmentation_data);
    return status;
}



extern "C" int seg_array_get_priors(PrecisionTYPE* prior_data)
{
    if (!NiftySeg_initialised())
        return 1;

    nifti_image * Prior = GlobSegmentation().SEG->GetPriorImage();
    if (Prior==NULL)
    {
         fprintf(stderr,"The prior image is not defined.\n");
         return 1;
    }

    int size_x = GlobSegmentation().SEG->get_size_x();
    int size_y = GlobSegmentation().SEG->get_size_y();
    int size_z = GlobSegmentation().SEG->get_size_z();  //fixme, only 3D
    int n_classes = GlobSegmentation().SEG->get_n_classes();

    memcpy((void*)prior_data, (void*)Prior->data, size_x*size_y*size_z*n_classes*sizeof(PrecisionTYPE));
     
    return 0;
}

extern "C" int seg_array_save_priors(char *filename)
{
    if (!NiftySeg_initialised())
         return 1;
    nifti_image *Prior=NULL;
    Prior = GlobSegmentation().SEG->GetPriorImage();
    if (Prior==NULL)
    {
         fprintf(stderr,"The prior image is not defined.\n");
         return 1;
    }
    nifti_set_filenames(Prior, filename, 0, 0);
    nifti_image_write(Prior);
    return 0;
}

extern "C" int seg_array_get_mask(PrecisionTYPE* mask_data)
{
    if (!NiftySeg_initialised())
        return 1;
    nifti_image *Mask=NULL;
    Mask = GlobSegmentation().SEG->GetMask();
    if (Mask==NULL)
    {
         fprintf(stderr,"The mask image is not defined.\n");
         return 1;
    }
    int size_x = GlobSegmentation().SEG->get_size_x();
    int size_y = GlobSegmentation().SEG->get_size_y();
    int size_z = GlobSegmentation().SEG->get_size_z();  //fixme, only 3D
    int n_classes = GlobSegmentation().SEG->get_n_classes();

    memcpy((void*)mask_data, (void*)Mask->data, size_x*size_y*size_z*n_classes*sizeof(bool));
    return 0;
}

extern "C" int seg_array_save_mask(char *filename)
{
    if (!NiftySeg_initialised())
         return 1;
    nifti_image *Mask=NULL;
    Mask = GlobSegmentation().SEG->GetMask();
    if (Mask==NULL)
    {
         fprintf(stderr,"The mask image is not defined.\n");
         return 1;
    }
    nifti_set_filenames(Mask, filename, 0, 0);
    nifti_image_write(Mask);
    return 0;
}


extern "C" int seg_array_get_Gaussian_parameters()
{
    if (!NiftySeg_initialised())
        return 1;

    return 0;
}


//extern "C" int seg_array_set_image(PrecisionTYPE *image, int size_x, int size_y, int size_z, int n_dimensions)
//{
//    nifti_image * T1;
//    seg_changeDatatype<PrecisionTYPE>(T1);
//}

extern "C" int seg_array_set_mask_fromfile(char *filename)
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().segment_param->filename_mask = filename;
    GlobSegmentation().segment_param->flag_mask=1;

    nifti_image * Mask=NULL;
    Mask = nifti_image_read(GlobSegmentation().segment_param->filename_mask,true);
    seg_changeDatatype<PrecisionTYPE>(Mask);
    if(Mask == NULL)
    {
        fprintf(stderr,"* Error when reading the mask image: %s\n",GlobSegmentation().segment_param->filename_mask);
        return 1;
    }
    GlobSegmentation().SEG->SetMaskImage(Mask);
    return 0; 
}


extern "C" int seg_array_set_priors_fromfiles(int n_classes, char **filenames)
{
    if (!NiftySeg_initialised())
        return 1;

    int status = 0;

    for(int classnum=0; classnum<GlobSegmentation().segment_param->numb_classes; classnum++)
        {
        GlobSegmentation().segment_param->filename_priors[classnum]=filenames[classnum];
        }
    GlobSegmentation().segment_param->flag_manual_priors=1;

    if (GlobSegmentation().segment_param->numb_classes != n_classes)
        {
        fprintf(stderr,"Inconsistent number of classes: init %d, input %d\n",GlobSegmentation().segment_param->numb_classes,n_classes);
        return 1;
        }
    nifti_image ** Priors_temp=new nifti_image * [GlobSegmentation().segment_param->numb_classes];
    nifti_image * Priors=NULL;

    for(int i=0; i<GlobSegmentation().segment_param->numb_classes; i++)
        {
        //fprintf(stderr,"NiftyRec: Reading Prior image: %s\n",GlobSegmentation().segment_param->filename_priors[i]);
        Priors_temp[i] = nifti_image_read(GlobSegmentation().segment_param->filename_priors[i],true);
        if(Priors_temp[i] == NULL)
            {
            fprintf(stderr,"* Error when reading the Prior image: %s\n",GlobSegmentation().segment_param->filename_priors[i]);
            return 1;
            }
        seg_changeDatatype<PrecisionTYPE>(Priors_temp[i]);
        }

    Priors=nifti_copy_nim_info(GlobSegmentation().SEG->GetInputImage());
    Priors->dim[0]=4;
    Priors->dim[4]=GlobSegmentation().segment_param->numb_classes;
    Priors->datatype=DT_FLOAT32;
    Priors->cal_max=1;

    nifti_update_dims_from_array(Priors);
    nifti_datatype_sizes(Priors->datatype,&Priors->nbyper,&Priors->swapsize);
    Priors->data = (void *) calloc(Priors->nvox, sizeof(PrecisionTYPE));

    Merge_Priors(Priors,Priors_temp,GlobSegmentation().segment_param);

    GlobSegmentation().SEG->SetPriorImage(Priors);

    for(int i=0;i<GlobSegmentation().segment_param->numb_classes;i++)
        {
        nifti_image_free(Priors_temp[i]);
        Priors_temp[i]=NULL;
        }
    delete [] Priors_temp;
GlobSegmentation().SEG->Init_EM();
    return status;
}

extern "C" int seg_array_get_loglikelihood(double *loglikelihood)
{
    if (!NiftySeg_initialised())
        return 1;
    *loglikelihood = GlobSegmentation().SEG->GetLoglik();
    return 0;
}

extern "C" int seg_array_get_mean_variance(double *mean, double *variance)
{
    if (!NiftySeg_initialised())
        return 1;

    int n_images = GlobSegmentation().SEG->get_n_images();
    int n_classes = GlobSegmentation().SEG->get_n_classes();
    float *mean_ptr = GlobSegmentation().SEG->GetMean();    
    float *variance_ptr = GlobSegmentation().SEG->GetVariance();    
    for (unsigned int i=0; i<n_images*n_classes;i++)
        mean[i]=mean_ptr[i];
    for (unsigned int i=0; i<n_images*n_images*n_classes;i++)
        variance[i]=variance_ptr[i];
    return 0;
}

extern "C" int seg_array_set_mean(float *mean)
{
    if (!NiftySeg_initialised())
        return 1;
    GlobSegmentation().SEG->SetMean(mean);
    return 0;
}

extern "C" int seg_array_set_variance(float *variance)
{
    if (!NiftySeg_initialised())
        return 1;
    GlobSegmentation().SEG->SetVariance(variance);
    return 0;
}


extern "C" int seg_array_cleanup()
{
    if (!NiftySeg_initialised())
        return 1;

    GlobSegmentation().initialised = false;
    
    GlobSegmentation().SEG->FreeInputImage();
    GlobSegmentation().SEG->FreeMaskImage();
    GlobSegmentation().SEG->FreeBiasImage();
    GlobSegmentation().SEG->FreePriorImage();

    return 0;
}




/* reg_maths functions - no persistent data */

extern "C" int seg_array_erode(float *image_ptr, int *image_size, float *out_ptr, int iterations)
{
     ImageSize * CurrSize = new ImageSize [1]();
     CurrSize->numel=image_size[0]*image_size[1]*image_size[2];
     CurrSize->xsize=image_size[0];
     CurrSize->ysize=image_size[1];
     CurrSize->zsize=image_size[2];
     CurrSize->usize=1;
     CurrSize->tsize=1;

     memcpy((void*) out_ptr, (void*) image_ptr, image_size[0]*image_size[1]*image_size[2]*sizeof(float));

     Erosion(out_ptr, iterations, CurrSize );
     delete [] CurrSize;
     return 0;
}



