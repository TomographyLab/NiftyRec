
function [image_out, lesions, radius, lesion_factor] = et_apply_lesions_2D(image_in, n_lesions, radius, lesion_factor)

im_size_x = size(image_in,1);
im_size_y = size(image_in,2);
image_out = image_in;
lesions = zeros(n_lesions,2);

for lesion = 1:n_lesions
    xc = rand*(im_size_x-2*radius)+radius;
    yc = rand*(im_size_y-2*radius)+radius;
    X = repmat((1:im_size_y),im_size_x,1);
    Y = repmat((1:im_size_x)',1,im_size_y);
    mask = find((X-xc).^2+(Y-yc).^2 < radius^2);
    image_out(mask) = image_out(mask)*lesion_factor;
    lesions(lesion,1) = xc;
    lesions(lesion,2) = yc;
end
