function wresult = savedcm(outputfilename,imaVOL,scaninfo)
% function wresult = saveecat(outputfilename,imaVOL,scaninfo)
%
% Matlab function to save ecat7 format output file. 
% The function use the ecatfile.m procedure coming from Flemming Hermansen.
%
% Matlab library function for MIA_gui utility. 
% University of Debrecen, PET Center/LB 2003
wresult = -1;

if scaninfo.Frames > 1
    hm = msgbox('Dynamic file cannot be saved.','MIA Info' );
    wresult = 0;
    return;
end

vol = permute(flipdim(flipdim(imaVOL,3),2),[2 1 3]);
%[ mh, hds, message ]   = ecatfile( 'blankhd', 'ecat7' );
hdin = load('ecat7_hdinfo.mat');
hd = hdin.hd;
hd.mh.patient_name= scaninfo.pnm;
hd.mh.patient_id = scaninfo.brn;
hd.mh.isotope_name = scaninfo.iso;
hd.mh.study_description = scaninfo.cntx;
hd.mh.file_system = 'ecat7';
hd.mh.file_type = 7;
hd.mh.num_planes = scaninfo.num_of_slice;
hd.mh.axial_fov = scaninfo.num_of_slice*scaninfo.pixsize(3)/10;
hd.mh.num_frames = 1;
hd.mh.num_gates= 1;
hd.mh.user_process_code = 'UDEB PETC';
[ fid, message ]       = ecatfile( 'create', outputfilename, hd.mh );
%message                = ecatfile( 'writemh', fid, mh );
%[ sh, hds, message ]   = ecatfile( 'blankhd', 'ecat7', 7 );
hd.sh.xyz_dimension(1:2)= scaninfo.imfm  ;
hd.sh.xyz_pixel_size= scaninfo.pixsize/10;
hd.sh.xyz_dimension(3) = scaninfo.num_of_slice ;
hd.sh.data_type = 6;
hd.sh.slice_width = scaninfo.pixsize(3)/10;
hd.sh.quant_scale = 1;
hd.sh.file_system = 'ecat7sh';
hd.sh.image_min = min(vol(:));
hd.sh.image_max = max(vol(:));
hd.sh.decay_corr_fctr = 1;
hd.sh.num_dimensions = length(size(vol));
selmatrix = [1 1 1 0 0];
message                = ecatfile( 'write', fid, vol, hd, selmatrix );
message                = ecatfile( 'close', fid );

wresult = 0;
