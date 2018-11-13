
% TT_DEMO_PROJECT
%     NiftyRec Demo: Transmission Tomography projection
%
%See also
%     TT_PROJECT, TT_BACKPROJECT
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


%% Load attenuation data 
attenuation = et_load_nifti('attenuation_01.nii'); 
attenuation = single(attenuation.img)*5e-7; 
attenuation(attenuation<=0)=0;

%% Parameters
detector_to_centre_mm        = 500;            % mm 
source_to_centre_rotation_mm = 1000;           % mm
N_projections                = 120;
detector_size_pix            = [180,480];      % pixels 
ray_step_vox                 = 0.3;            % voxels
volume_size_mm               = [416,318,181];  % mm
detector_size_mm_0           = 1.5*[180,480];  % mm 
rotation_angle_rad           = 2*pi;           % rad

%% Define scanner geometry 
% Initial geometry
center_rotation_mm        = [volume_size_mm(1)/2,volume_size_mm(2)/2,volume_size_mm(3)/2];
source_location_mm_0      = [source_to_centre_rotation_mm+center_rotation_mm(1), center_rotation_mm(2), center_rotation_mm(3)]; % mm
detector_rotation_rad_0   = [0, -pi/2, 0];  % radians
detector_translation_mm_0 = [center_rotation_mm(1)-detector_to_centre_mm, (center_rotation_mm(1)-detector_size_mm_0(2))/2, (center_rotation_mm(2)-detector_size_mm_0(1))/2]; 

detector_size_mm        = ones(N_projections,2);
source_location_mm      = zeros(N_projections,3);
detector_translation_mm = zeros(N_projections,3);
detector_rotation_rad   = zeros(N_projections,3);

% Rotate the detector and the source around the object 
rotation_rad = 0;
for i=1:N_projections
    detector_size_mm(i,:)  = detector_size_mm_0;
    detector_rotation_rad(i,:) = detector_rotation_rad_0;

    rotation_rad = rotation_rad + rotation_angle_rad / N_projections; 
    detector_rotation_rad(i,3) = rotation_rad;
    M = [cos(rotation_rad), -sin(rotation_rad), 0; sin(rotation_rad), cos(rotation_rad), 0; 0,0,1];
    detector_translation_mm(i,:)  = (M*(detector_translation_mm_0'-center_rotation_mm'))'+center_rotation_mm;
    source_location_mm(i,:) = (M*(source_location_mm_0'-center_rotation_mm'))'+center_rotation_mm;
end

p = tt_project_ray_mex(single(attenuation),volume_size_mm,source_location_mm,int32(detector_size_pix),detector_size_mm,detector_translation_mm,detector_rotation_rad,ray_step_vox);
p = exp(-p); 
min_invp = min(1./p(:)); 
max_invp = max(1./p(:)); 

figure; image((-min_invp+1./p(:,:,1))/max_invp*128); colormap gray; axis tight off equal; 
for i=1:N_projections
    image((-min_invp+1./p(:,:,i))/max_invp*128); colormap gray; axis equal tight off; pause(0.1);
end

