/*
 *  rec_spect.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <tclap/CmdLine.h>
#include "niftyrec_version.h"
#include "_et_array_interface.h"

#ifdef _SUPPORT_NRRD
#include "teem/nrrd.h"
#endif 

#define FILETYPE_NIFTI 0
#define FILETYPE_NRRD  1

#define BACKGROUND_ACTIVITY 0.0f
#define BACKGROUND_ATTENUATION 0.0f
#define EPSILON 0.0000001f

#define XML_FILE "./rec_spect.xml"

char xml_string[] = ""
"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
"<executable>\n"
"  <category>Reconstruction</category>\n"
"  <title>NiftyRec SPECT</title>\n"
"  <description>\n"
"  Reconstruct SPECT images from sinogram data\n"
"  </description>\n"
"  <version>1.5</version>\n"
"  <documentation-url>http://sourceforge.net/niftyrec</documentation-url>\n"
"  <license>http://sourceforge.net/niftyrec</license>\n"
"  <contributor>CMIC-UCL</contributor>\n"
"  <parameters>\n"
"    <label>Workspace IO</label>\n"
"    <description>Input images</description>\n"

"    <image>\n"
"      <name>Sinogram</name>\n"
"      <label>Sinogram</label>\n"
"      <channel>input</channel>\n"
"      <longflag>sinogram</longflag>\n"
"      <description>Sinogram data</description>\n"
"    </image>\n"

"    <image>\n"
"      <name>Attenuation</name>\n"
"      <label>Attenuation</label>\n"
"      <channel>input</channel>\n"
"      <longflag>attenuation</longflag>\n"
"      <description>Attenuation map</description>\n"
"    </image>\n"

"    <image>\n"
"      <name>OutputActivity</name>\n"
"      <label>Output Activity</label>\n"
"      <channel>output</channel>\n"
"      <longflag>output</longflag>\n"
"      <description>Output activity</description>\n"
"    </image>\n"

"  </parameters>\n"
"  <parameters>\n"
"    <label>SPECT Acquisition Parameters</label>\n"
"    <description>Parameters of the SPECT acquisition</description>\n"

"    <double>\n"
"      <name>FirstCamera</name>\n"
"      <longflag>firstcamera</longflag>\n"
"      <description>\n"
"      First camera position in degrees\n"
"      </description>\n"
"      <label>First Camera Position</label>\n"
"      <default>0.0</default>\n"
"    </double>\n"

"    <double>\n"
"      <name>LastCamera</name>\n"
"      <longflag>lastcamera</longflag>\n"
"      <description>\n"
"      Last camera position in degrees\n"
"      </description>\n"
"      <label>Last Camera Position</label>\n"
"      <default>180.0</default>\n"
"    </double>\n"

"    <integer>\n"
"      <name>RotationAxis</name>\n"
"      <longflag>axis</longflag>\n"
"      <description>Axis of rotation of the Gamma Camera</description>\n"
"      <label>RotationAxis</label>\n"
"      <default>0</default>\n"
"      <constraints>\n"
"        <minimum>0</minimum>\n"
"        <maximum>2</maximum>\n"
"        <step>1</step>\n"
"      </constraints>\n"
"    </integer>\n"

"    <double>\n"
"      <name>FWHM0</name>\n"
"      <longflag>fwhm0</longflag>\n"
"      <description>\n"
"      PSF Full Width Half Maximum in pixels at distance dist0\n"
"      </description>\n"
"      <label>PSF Full Width Half Maximum 0</label>\n"
"      <default>3.0</default>\n"
"      <constraints>\n"
"        <minimum>0.0</minimum>\n"
"        <maximum>256.0</maximum>\n"
"      </constraints>\n"
"    </double>\n"

"    <double>\n"
"      <name>FWHM1</name>\n"
"      <longflag>fwhm1</longflag>\n"
"      <description>\n"
"      PSF Full Width Half Maximum in pixels at distance dist1\n"
"      </description>\n"
"      <label>PSF Full Width Half Maximum 1</label>\n"
"      <default>3.0</default>\n"
"      <constraints>\n"
"        <minimum>0.0</minimum>\n"
"        <maximum>256.0</maximum>\n"
"        <step>0.1</step>\n"
"      </constraints>\n"
"    </double>\n"

"    <double>\n"
"      <name>Efficiency0</name>\n"
"      <longflag>efficiency0</longflag>\n"
"      <description>\n"
"      PSF efficiency at distance dist0\n"
"      </description>\n"
"      <label>PSF Efficiency 0</label>\n"
"      <default>0.9</default>\n"
"      <constraints>\n"
"        <minimum>0.0</minimum>\n"
"        <maximum>1.0</maximum>\n"
"        <step>0.05</step>\n"
"      </constraints>\n"
"    </double>\n"

"    <double>\n"
"      <name>Efficiency1</name>\n"
"      <longflag>efficiency1</longflag>\n"
"      <description>\n"
"      PSF efficiency at distance dist1\n"
"      </description>\n"
"      <label>PSF Efficiency 1</label>\n"
"      <default>0.9</default>\n"
"      <constraints>\n"
"        <minimum>0.0</minimum>\n"
"        <maximum>1.0</maximum>\n"
"        <step>0.05</step>\n"
"      </constraints>\n"
"    </double>\n"

"    <double>\n"
"      <name>Distance0</name>\n"
"      <longflag>dist0</longflag>\n"
"      <description>\n"
"      Distance from the camera in pixels of PSF0 (characterised by FWHM0 and Efficiency0)\n"
"      </description>\n"
"      <label>PSF distance 0</label>\n"
"      <default>10.0</default>\n"
"      <constraints>\n"
"        <minimum>0.0</minimum>\n"
"        <maximum>4096.0</maximum>\n"
"        <step>0.1</step>\n"
"      </constraints>\n"
"    </double>\n"

"    <double>\n"
"      <name>Distance1</name>\n"
"      <longflag>dist1</longflag>\n"
"      <description>\n"
"      Distance from the camera in pixels of PSF1 (characterised by FWHM1 and Efficiency1)\n"
"      </description>\n"
"      <label>PSF distance 1</label>\n"
"      <default>10.0</default>\n"
"      <constraints>\n"
"        <minimum>0.0</minimum>\n"
"        <maximum>4096.0</maximum>\n"
"        <step>0.1</step>\n"
"      </constraints>\n"
"    </double>\n"

"    <double>\n"
"      <name>RotationRadius</name>\n"
"      <longflag>radius</longflag>\n"
"      <description>\n"
"      Radius of rotation of the Gamma Camera in pixels\n"
"      </description>\n"
"      <label>Rotation Radius</label>\n"
"      <default>20.0</default>\n"
"      <constraints>\n"
"        <minimum>0.0</minimum>\n"
"        <maximum>4096.0</maximum>\n"
"        <step>0.1</step>\n"
"      </constraints>\n"
"    </double>\n"

"  </parameters>\n"
"  <parameters>\n"
"    <label>Reconstruction Parameters</label>\n"
"    <description>\n"
"    Parameters for the reconstruction\n"
"    </description>\n"

"    <integer>\n"
"      <name>iterations</name>\n"
"      <longflag>iterations</longflag>\n"
"      <description>\n"
"      Iterations of the reconstruction algorithm.\n"
"      </description>\n"
"      <label>Iterations</label>\n"
"      <default>20</default>\n"
"      <constraints>\n"
"        <minimum>1</minimum>\n"
"        <maximum>5000</maximum>\n"
"        <step>1</step>\n"
"      </constraints>\n"
"    </integer>\n"

"    <boolean>\n"
"      <name>ForceGpuAccelerationOFF</name>\n"
"      <longflag>gpu_off</longflag>\n"
"      <description>\n"
"      Force OFF GPU acceleration. Otherwise NiftyRec will attempt to use GPU hardware acceleration. \n"
"      </description>\n"
"      <label>Force disable GPU</label>\n"
"      <default>false</default>\n"
"    </boolean>\n" 

"    <integer>\n"
"      <name>SubsetsOSEM</name>\n"
"      <longflag>subsets</longflag>\n"
"      <description>\n"
"      Number of subsets for OSEM reconstruction. The higher, the faster, the lower, the more accurate. \n"
"      </description>\n"
"      <label>Subsets OSEM</label>\n"
"      <default>16</default>\n"
"      <constraints>\n"
"        <minimum>1</minimum>\n"
"        <maximum>180</maximum>\n"
"        <step>1</step>\n"
"      </constraints>\n"
"    </integer>\n"

"  </parameters>\n"
"</executable>\n";




int stream_xml()
{
        //if the XML file is in the path, read it and stream it to stdout, otherwise use default from compile time
        std::string line;
        std::ifstream xml_file (XML_FILE);
        if (xml_file.is_open())
        {
            while ( xml_file.good() )
            {
                std::getline (xml_file,line);
                std::cout << line << std::endl;
            }
        xml_file.close();
        }
        else 
        {
        std::cout << xml_string;
        }
        return 0;   
}

#ifdef _SUPPORT_NRRD
nifti_image *nrrd_to_nifti(Nrrd *nrrdImage)
{
    nifti_image *niftiImage = NULL;
    /* say something about the array */
    std::cout << "================ Nrrd -> Nifti conversion ===============" << std::endl;
    printf("-- Converting %d-dimensional NRRD of type %d (%s)\n", nrrdImage->dim, nrrdImage->type, airEnumStr(nrrdType, nrrdImage->type));
    for (int n_axis=0; n_axis<nrrdImage->dim; n_axis++)
        printf("-- Size axis %d: %d\n",n_axis,nrrdImage->axis[n_axis].size);
    printf("-- The array contains %d elements, each %d bytes in size\n", (int)nrrdElementNumber(nrrdImage), (int)nrrdElementSize(nrrdImage));
    /* convert it: size and dimensions */
    int dim[8];
    dim[0] = nrrdImage->dim;
    for (int n_axis=0; n_axis<nrrdImage->dim; n_axis++)
        dim[n_axis+1]=nrrdImage->axis[n_axis].size;
    for (int n_axis=nrrdImage->dim+1; n_axis<8; n_axis++)
        dim[n_axis]=1;
    /* convert it: spatial transformations */
    //.. set axis.measurementFrame, axis.spaceOrigin, axis.spaceUnits
    if (nrrdImage->type == nrrdTypeFloat) 
    {
        //allocate new Nifti image struct
        niftiImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        //copy data (rather than passing the pointer, to simplify free - use nrrdNuke() and nifti_image_free() )
        niftiImage->data = (void *) calloc((int)nrrdElementNumber(nrrdImage),sizeof(float));
        memcpy((void*)niftiImage->data, (void*)nrrdImage->data,(int)nrrdElementNumber(nrrdImage)*sizeof(float));
        std::cout << "-- Conversion done" << std::endl;
    }
    else
    {
        std::cout << "-- Conversion from Nrrd to Nifti not supported for Nrrd format %s." << airEnumStr(nrrdType, nrrdImage->type) << std::endl;
    }
    std::cout << "=========================================================" << std::endl;
    return niftiImage;
}


Nrrd *nifti_to_nrrd(nifti_image *niftiImage)
{
    std::cout << "================ Nifti -> Nrrd conversion ===============" << std::endl;
    Nrrd *nrrdImage = nrrdNew(); 
    unsigned int dim = niftiImage->dim[0];
    size_t size[NRRD_DIM_MAX];
    for (int n_axis=0; n_axis<dim; n_axis++)
        size[n_axis] = niftiImage->dim[n_axis+1];
    if(niftiImage->datatype == NIFTI_TYPE_FLOAT32)
    {
        if (nrrdAlloc_nva(nrrdImage, nrrdTypeFloat, dim, size)) 
        {
            std::cout << "-- Conversion failed, cannot allocate the NRRD image. " << std::endl;
            return NULL;
        } 
        /* copy data - rather than passing the pointer, to facilitate a bit memory handling */
        memcpy((void*)nrrdImage->data,(void*)niftiImage->data,(int)nrrdElementNumber(nrrdImage)*sizeof(float));
    }
    else
    {
        std::cout << "-- Conversion from Nifti to Nrrd not supported for Nifti data type %d." << niftiImage->datatype << std::endl;
        return NULL;
    }
    printf("-- Created %d-dimensional NRRD of type %d (%s)\n", nrrdImage->dim, nrrdImage->type, airEnumStr(nrrdType, nrrdImage->type));
    for (int n_axis=0; n_axis<nrrdImage->dim; n_axis++)
        printf("-- Size axis %d: %d\n",n_axis,nrrdImage->axis[n_axis].size);
    std::cout << "-- Conversion done" << std::endl;
    std::cout << "=========================================================" << std::endl;
    return nrrdImage;
}
#endif




int main(int argc, char** argv)
{
    int status;
    if (argc==2)
    {
        if(strcmp(argv[1], "--xml") == 0)
        {
            stream_xml();
            return 0;
        }
    }       
        
    try 
    {  
    //Parse command line arguments
        TCLAP::CmdLine cmd("NiftyRec SPECT - Maximum Likelihood Expectation Maximisation (MLEM) and Odered Subsets Expectation Maximisation (OSEM) reconstruction of SPECT sinogram data. ", ' ', VERSION);

        //Required input sinogram data file
        TCLAP::ValueArg<std::string> sinogramArg("s","sinogram","Sinogram data file to reconstruct. Volumetric data of size [NxNxn_cameras].",true,"none","string",cmd);
        //Required output file name
        TCLAP::ValueArg<std::string> outputArg("o","output","File name for the output reconstructed image. Reconstructed image will be of size [NxMxN].",true,"none","string",cmd);
        //Optional attenuation map
        TCLAP::ValueArg<std::string> attenuationArg("a","attenuation","Filename of the attenuation map image. Image must be of size [NxMxN] for sinogram of size [NxN]. ",false,"none","string",cmd);
        //Optional first camera angle
        TCLAP::ValueArg<float> firstcameraArg("c","firstcamera","Position of the first camera in degrees. Defaults to 0.0 deg",false,0.0f,"float",cmd);
        //Optional last camera angle
        TCLAP::ValueArg<float> lastcameraArg("C","lastcamera","Position of the last camera in degrees. Defaults to 180.0 deg",false,180.0f,"float",cmd);
        //Optional axis of rotation
        TCLAP::ValueArg<int> axisArg("x","axis","Axis of rotation of the Gamma Camera [0,1,2]. Defaults to 0.",false,0,"int",cmd);
        //Optional FWHM0 
        TCLAP::ValueArg<float> fwhm0Arg("f","fwhm0","FWHM of the PSF at distance dist0 in pixels",false,3.0f,"float",cmd);
        //Optional FWHM1 
        TCLAP::ValueArg<float> fwhm1Arg("F","fwhm1","FWHM of the PSF at distance dist1 in pixels",false,3.0f,"float",cmd);
        //Optional efficiency0
        TCLAP::ValueArg<float> efficiency0Arg("e","efficiency0","Efficiency of the collimator/detector response (PSF) at distance dist0",false,0.9f,"float",cmd);
        //Optional efficiency1
        TCLAP::ValueArg<float> efficiency1Arg("E","efficiency1","Efficiency of the collimator/detector response (PSF) at distance dist1",false,0.9f,"float",cmd); 
        //Optional dist0
        TCLAP::ValueArg<float> dist0Arg("d","dist0","Distance of the PSF characterisation 0 in pixels (also specify fwhm0 and efficiency0)",false,10.0f,"float",cmd); 
        //Optional dist1
        TCLAP::ValueArg<float> dist1Arg("D","dist1","Distance of the PSF characterisation 1 in pixels (also specify fwhm1 and efficiency1)",false,10.0f,"float",cmd); 
        //Optional rotation radius
        TCLAP::ValueArg<float> radiusArg("r","radius","Radius of rotation of the Gamma camera in pixels. Defaults to 20.0",false,20.0f,"float",cmd); 
        //Optional number of iterations
        TCLAP::ValueArg<int> iterationsArg("i","iterations","MLEM iterations. Defaults to 20",false,20,"int",cmd);
        //Optional number of OSEM subsets
        TCLAP::ValueArg<int> subsetsArg("u","subsets","Number of subsets for OSEM reconstruction. Defaults to 16. Set to 0 or 1 (indifferently) for MLEM. ",false,16,"int",cmd);
        //Optional switch use gpu
        TCLAP::SwitchArg gpuSwitch("g","gpu_off","[1] - Graphics Processing Unit (GPU) hardware acceleration OFF. [0] - GPU ON. Defaults to 0.", cmd, false);
        //Optional switch verbose
        TCLAP::SwitchArg verboseSwitch("v","verbose","Print verbose information", cmd, false);

        cmd.parse( argc, argv );

        std::string sinogram_filename = sinogramArg.getValue().c_str();
        std::string output_filename = outputArg.getValue().c_str();
        std::string attenuation_filename = attenuationArg.getValue().c_str();
        float firstcamera = firstcameraArg.getValue(); 
        float lastcamera = lastcameraArg.getValue(); 
        unsigned int axis = axisArg.getValue(); 
        float PSF_fwhm0 = fwhm0Arg.getValue(); 
        float PSF_fwhm1 = fwhm1Arg.getValue(); 
        float PSF_dist0 = dist0Arg.getValue(); 
        float PSF_dist1 = dist1Arg.getValue(); 
        float PSF_efficiency0 = efficiency0Arg.getValue(); 
        float PSF_efficiency1 = efficiency1Arg.getValue(); 
        float radius = radiusArg.getValue(); 
        int iterations = iterationsArg.getValue();
        int gpu = 1; if (gpuSwitch.getValue()) gpu=0; 
        int subsets = subsetsArg.getValue(); 
        bool verbose = verboseSwitch.getValue(); 

        //Load input sinogram
        int filetype_sinogram;
        nifti_image *sinogramImage = NULL;

        if (sinogram_filename.substr(sinogram_filename.find_last_of(".") + 1) == "nii" || sinogram_filename.substr(sinogram_filename.find_last_of(".") + 1) == "nii.gz")
        {
            std::cout << "Nifti sinogram file: " <<sinogram_filename<< std::endl;
            filetype_sinogram = FILETYPE_NIFTI;
        }
        else if (sinogram_filename.substr(sinogram_filename.find_last_of(".") + 1) == "nrrd")
        {
            std::cout << "NRRD sinogram file: " <<sinogram_filename<< std::endl;
            filetype_sinogram = FILETYPE_NRRD;
        }
        else
        {
            std::cout << "Unknown file format" << std::endl;
            return 1;
        }

        if (filetype_sinogram == FILETYPE_NIFTI)
        {
            sinogramImage = nifti_image_read(sinogram_filename.c_str(),true);
            if (sinogramImage == NULL)
            {
                std::cout << "Couldn't load the sinogram " << sinogram_filename << std::endl;
                return 1;                
            }
        }
        else if (filetype_sinogram == FILETYPE_NRRD)
        {
#ifdef _SUPPORT_NRRD
            /* create a nrrd; at this point this is just an empty container */
            Nrrd *sinogram_nrrdImage;
            sinogram_nrrdImage = nrrdNew();
            char *err;

            /* read in the nrrd from file */
            if (nrrdLoad(sinogram_nrrdImage, sinogram_filename.c_str(), NULL))
            {
                err = biffGetDone(NRRD);
                fprintf(stderr, "Problem reading \"%s\":\n%s", sinogram_filename.c_str(), err);
                free(err);
                return 1;
            }
            /* convert NRRD to Nifti */

            sinogramImage = nrrd_to_nifti(sinogram_nrrdImage);
            if (sinogramImage == NULL)
            {
                std::cout << "Convertion from Nrrd to Nifti failed." << std::endl;
                return 1;
            }

            /* blow away both the Nrrd struct *and* the memory at nin->data. (nrrdNix() frees the struct but not the data, nrrdEmpty() frees the data but not the struct) */
            nrrdNuke(sinogram_nrrdImage);
#else 
            std::cout << "NRRD file format not supported. Please enable the NRRD support flag at compile time." << std::endl;
            return 1;
#endif 
        }
        else
        {
            std::cout << "Unknown file format" << std::endl;
            return 1;            
        }


        //Reconstruct
        float *sinogram_data = (float*)sinogramImage->data;
        unsigned int size_x = sinogramImage->nx;
        unsigned int size_y = sinogramImage->ny;
        int n_cameras = sinogramImage->nz;
        unsigned int psf_size_x; 
        unsigned int psf_size_y;
        const int use_psf = 1;
        int use_attenuation = 0; 
        float *activity_data = (float*) malloc(size_x*size_y*size_x*sizeof(float)); 
        float *psf_data=NULL; 
        float *attenuation_data=NULL; 
        if (1)
        {
            std::cout << "=============== NiftyRec SPECT Parameters ===============" << std::endl;
            std::cout << "Sinogram file name:       " << sinogram_filename << std::endl;
            std::cout << "Output file name:         " << output_filename << std::endl;
            std::cout << "Attenuation file name:    " << attenuation_filename << std::endl;
            std::cout << "First camera:             " << firstcamera << std::endl;
            std::cout << "Last camera:              " << lastcamera << std::endl;
            std::cout << "Rotation axis:            " << axis << std::endl;
            std::cout << "PSF distance 0 [pixels]:  " << PSF_dist0 << std::endl;
            std::cout << "PSF FWHM 0 [pixels]:      " << PSF_fwhm0 << std::endl;
            std::cout << "PSF efficiency 0:         " << PSF_efficiency0 << std::endl;
            std::cout << "PSF distance 1 [pixels]:  " << PSF_dist1 << std::endl;
            std::cout << "PSF FWHM 1 [pixels]:      " << PSF_fwhm1 << std::endl;
            std::cout << "PSF efficiency 1:         " << PSF_efficiency1 << std::endl;
            std::cout << "Rotation radius [pixels]: " << radius << std::endl;
            std::cout << "Iterations:               " << iterations << std::endl;
            std::cout << "OSEM subsets:             " << subsets << std::endl;
            std::cout << "Use GPU:                  " << gpu << std::endl;
            std::cout << "Sinogram size X [pixels]: " << size_x << std::endl; 
            std::cout << "Sinogram size Y [pixels]: " << size_y << std::endl; 
            std::cout << "Number of projections:    " << n_cameras << std::endl; 
            std::cout << "Activity size X [pixels]: " << size_x << std::endl; 
            std::cout << "Activity size Y [pixels]: " << size_y << std::endl; 
            std::cout << "Activity size Z [pixels]: " << size_x << std::endl;
            std::cout << "=========================================================" << std::endl;
        }

        status = et_array_calculate_size_psf(&psf_size_x, &psf_size_y, PSF_fwhm0, PSF_efficiency0, PSF_dist0, PSF_fwhm1, PSF_efficiency1, PSF_dist1); 
        if (status != 0) { fprintf(stderr,"rec_spect: Error creating PSF.\n"); return status; }
        psf_data = (float*) malloc(psf_size_x*psf_size_y*size_x*sizeof(float)); 
        status = et_array_make_psf(psf_data, psf_size_x, psf_size_y, PSF_fwhm0, PSF_efficiency0, PSF_dist0, PSF_fwhm1, PSF_efficiency1, PSF_dist1, size_x); 
        if (status != 0) { fprintf(stderr,"rec_spect: Error creating PSF.\n"); free(psf_data); return status; }

        if (use_attenuation) 
            {
            attenuation_data = (float*) malloc(size_x*size_y*size_x*sizeof(float)); 
            } 

        if (subsets<=1)
            status = et_array_mlem(activity_data, size_x, size_y, sinogram_data, n_cameras, firstcamera, lastcamera, axis, iterations, use_attenuation, attenuation_data, use_psf, psf_data, psf_size_x, psf_size_y, BACKGROUND_ACTIVITY, BACKGROUND_ATTENUATION, EPSILON, gpu);
        else
            status = et_array_osem(activity_data, size_x, size_y, subsets, sinogram_data, n_cameras, firstcamera, lastcamera, axis, iterations, use_attenuation, attenuation_data, use_psf, psf_data, psf_size_x, psf_size_y, BACKGROUND_ACTIVITY, BACKGROUND_ATTENUATION, EPSILON, gpu);

        if (status!=0) 
            { 
            fprintf(stderr,"rec_spect: reconstruction failed.\n "); 
            free(activity_data); 
            if (use_psf)
                free(psf_data);
            if (use_attenuation)
                free(attenuation_data); 
            return 1; 
            } 

        //Save result
        int filetype_output;

        if (output_filename.substr(output_filename.find_last_of(".") + 1) == "nii" || output_filename.substr(output_filename.find_last_of(".") + 1) == "nii.gz")
        {
            std::cout << "Nifti output file: " <<output_filename<< std::endl;
            filetype_output = FILETYPE_NIFTI;
        }
        else if (output_filename.substr(output_filename.find_last_of(".") + 1) == "nrrd")
        {
            std::cout << "NRRD output file: " <<output_filename<< std::endl;
            filetype_output = FILETYPE_NRRD;
        }
        else
        {
            std::cout << "Unknown file format" << std::endl;
            return 1;
        }

        int dim[8];
        dim[0]=3;
        dim[1]=size_x;
        dim[2]=size_y;
        dim[3]=size_x;
        dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1;
        nifti_image *activityImage;
        activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        activityImage->data = static_cast<void *>(activity_data);        
        nifti_set_filenames(activityImage, output_filename.c_str(), 0, 0);

        if (filetype_output == FILETYPE_NIFTI)
        {
            nifti_image_write(activityImage);
        }
        else if (filetype_output == FILETYPE_NRRD)
        {
#ifdef _SUPPORT_NRRD
            /* write out the nrrd to a different file */
            Nrrd *activity_nrrdImage = nifti_to_nrrd(activityImage);
            char *err;
//            nrrdLoad(activity_nrrdImage, sinogram_filename.c_str(), NULL);//
            if (nrrdSave(output_filename.c_str(), activity_nrrdImage, NULL))
            {
                err = biffGetDone(NRRD);
                fprintf(stderr, "Trouble writing \"%s\":\n%s", output_filename.c_str(), err);
                free(err);
                return 1;
            }
            nrrdNuke(activity_nrrdImage);
#else
            std::cout << "NRRD file format not supported. Please enable the NRRD support flag at compile time." << std::endl;
            return 1;
#endif
        }
        //Clean up
        nifti_image_free(activityImage);

    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    } 
    return 0;
}


