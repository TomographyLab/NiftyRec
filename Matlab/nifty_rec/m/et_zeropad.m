function out_image = et_zeropad(in_image,N_x,N_y,N_z)
   [x,y,z] = size(in_image);
   im(floor((max_size-x)/2)+1:floor((max_size-x)/2)+x,floor((max_size-y)/2)+1:floor((max_size-y)/2)+y,floor((max_size-z)/2)+1:floor((max_size-z)/2)+z) = mri; 
end