function [imaVOL, scaninfo, hd] = loadecat(filename)
%function [imaVOL,scaninfo,hd] = loadecat(filename)
% This function loads ecat6.4 or ecat7(even dynamic) format input files. 
% On reading it uses the ecatfile.m, readecatvol.m procedures coming 
% from Flemming Hermansen (Aarhus University Hospitals).
%
% Matlab library function for MIA_gui utility. 
% University of Debrecen, PET Center/LB 2003


if nargin == 0
     [FileName, FilePath] = uigetfile('*.v','Select ecat file');
     filename = [FilePath,FileName];
     if FileName == 0;
          imaVOL = [];scaninfo = [];
          return;
     end
end

ecat_datain = readecatvol( filename, [ -inf, -inf, -inf, -inf, -inf; inf, inf, inf, inf, inf ] );
%[ fid, message1 ]       = ecatfile( 'open', filename );
%[ matranges, message2 ] = ecatfile( 'matranges', fid );
%[ vol, hd, message3 ]   = ecatfile( 'read', fid, matranges(1,:)  );
%message4                = ecatfile( 'close', fid );
hd = ecat_datain.hd{1};


% get info from minc file and save them in the scaninfo structure 
%
scaninfo.pnm	   = hd.mh.patient_name;
scaninfo.brn       = hd.mh.patient_id;
scaninfo.rid	   = [];
scaninfo.rin	   = [];	
scaninfo.daty	   = [];
scaninfo.datm	   = [];
scaninfo.datd      = [];
scaninfo.timh	   = [];
scaninfo.timm	   = [];
scaninfo.tims	   = [];
scaninfo.mtm       = [];
scaninfo.iso 	   = hd.mh.isotope_name;
scaninfo.half      = [];
scaninfo.trat      = [];
%scaninfo.imfm  	   = hd.sh.xyz_dimension(1:2);
scaninfo.cntx      = hd.mh.study_description;
scaninfo.cal       = [];
scaninfo.min       = [];
scaninfo.mag       = [];
%scaninfo.pixsize =  hd.sh.xyz_pixel_size*10;%mm
scaninfo.start_times     =  ecat_datain.times(:,1)';
scaninfo.frame_lengths    =  ecat_datain.times(:,4)';
scaninfo.tissue_ts  =  ecat_datain.times(:,2)';
scaninfo.Frames     = size(ecat_datain.times(:,1),1);
%scaninfo.num_of_slice    = hd.sh.xyz_dimension(3);
scaninfo.FileType    = hd.mh.file_system;
scaninfo.float = 0;

if  scaninfo.Frames == 1
    vol = ecat_datain.vol{1};
    imaVOL = flipdim(flipdim(permute(vol,[2 1 3]),2),3);
else
    imaVOL = int16(zeros(scaninfo.imfm(1),scaninfo.imfm(2),scaninfo.num_of_slice,scaninfo.Frames));
    for i=1:scaninfo.Frames
        vol = ecat_datain.vol{i};
        imaVOL(:,:,:,i) =  flipdim(flipdim(permute(vol,[2 1 3]),2),3);
    end
end

% zero padding the negativ elements
imaVOL(find(imaVOL < 0)) = 0;

