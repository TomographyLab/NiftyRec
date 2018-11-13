function [joint_entropy, varargout] = et_joint_entropy_nonparametric(input_data_1, input_data_2, varargin)
%ET_JOINT_ENTROPY_NONPARAMETRIC 
%    Fast and accurate bivariate kernel density estimator with Gaussian kernel and diagonal
%    bandwidth matrix. The function returns the Joint Entropy between the two variables measured 
%    on the kernel density estimate; the gradient of the Joint Entropy with
%    respect of the input variables and, optionally, the bivariate density estimate. 
%    The function also optionally returns the marginal densities, the
%    marginal entropies and the gradient of the marginal entropies with
%    respect to each of the input variables. 
%
%Description
%
%    [JOINT_ENTROPY, GRADIENT_JOINT_ENTROPY_1, JOINT_DENSITY, JOINT_HISTOGRAM, 
%    HIST_LIMITS_1, HIST_LIMITS_2, GRADIENT_JOINT_ENTROPY_2, ENTROPY_1, ENTROPY_2, 
%    GRADIENT_ENTROPY_1, GRADIENT_ENTROPY_2] = et_jointentropy_nonparametric(INPUT_DATA_1, 
%    INPUT_DATA_2, N_BINS_1, N_BINS_2, BANDWIDTH_1, BANDWIDTH_2,
%    HIST_LIMITS_1, HIST_LIMITS_2); 
%
%  All input arguments are optional except for the input data INPUT_DATA_1
%  and INPUT_DATA_2. 
%  All output arguments are optional except for the first (JOINT ENTROPY). 
%
%
%  Input arguments: 
%
%    INPUT_DATA_1: this is a [Nx1] array of continuous values. 
%       INPUT_DATA_1 can also be specified as a multi-dimensional array (e.g. 2D image,
%       3D image). In this case the gradient (output) will
%       have the same size as the input data. The size of INPUT_DATA_1 must 
%       be equal to the size of INPUT_DATA_2. 
%
%    INPUT_DATA_2: this should be a [Nx1] array of continuous values.
%       INPUT_DATA_2 can also be specified as a multi-dimensional array (e.g. 2D image,
%       3D image). In this case the gradient, returned by the function, will
%       have the same size as the input data. The size of INPUT_DATA_2 must 
%       be equal to the sizeof INPUT_DATA_1. 
%
%
%  Optional input arguments: 
%
%    N_BINS_1, N_BINS_2: number of bins for the discretisation of the joint
%       histogram, along the axis INPUT_DATA_1 and INPUT_DATA_2. 
%       The joint histogram is constructed on a square grid of
%       N_BINS_1*N_BINS_2 elements. 
%       Construction of the joint histogram and smoothing are a computationally
%       efficient way to discretise the kernel density estimation problem. 
%       The higher the value of N_BINS_1 and N_BINS_2, the more accurate the estimate of the
%       joint density. The lower the value of N_BINS_1 and N_BINS_2, the faster the
%       algorithm and the lower the memory requirements. 
%       Default value: 256. Typical values: 256 or 512. 
%
%    BANDWIDTH_1, BANDWIDTH_2: bandwidths (Full Width at Half Maximum - FWHM) of the 
%       Gaussian kernel for variables INPUT_DATA_1 and INPUT_DATA_2. 
%       FWHM=1 means that the FWHM of the kernel is as large as the domain of
%       the histogram. 
%       Typical value: 0.02 - 0.2 
%       If not set, then the bandwidth is selected automatically. 
%
%    HIST_LIMITS_1: [min_1,max_1] minimum and maximum values of variable 
%       INPUT_DATA_1 for the construction of the joint histogram. 
%       They are autimatically set if not specified. Values that fall outside 
%       of the range are histogrammed in the first and last bins. 
%
%    HIST_LIMITS_2: [min_2,max_2] minimum and maximum values of variable 
%       INPUT_DATA_2 for the construction of the joint histogram. 
%       They are autimatically set if not specified. Values that fall outside 
%       of the range are histogrammed in the first and last bins.  
%
%
%  Output arguments: 
%
%    JOINT_ENTROPY Joint Entropy of INPUT_DATA_1 and INPUT_DATA_2, 
%       calculated according to the kernel joint density estimate. 
%       It is a scalar value. 
% 
%
%  Optional output arguments: 
%
%    GRADIENT_JOINT_ENTROPY_1 gradient of the Joint Entropy with respect 
%       to INPUT_DATA_1. GRADIENT_JOINT_ENTROPY_1 has the same size as
%       INPUT_DATA_1 (and INPUT_DATA_2). 
%
%    JOINT_DENSITY: kernel density estimate of (INPUT_DATA_1,INPUT_DATA_2).
%
%    JOINT_HISTOGRAM: joint histogram of (INPUT_DATA_1,INPUT_DATA_2). 
% 
%    HIST_LIMITS_1: [min_1,max_1] minimum and maximum values of variable 
%       INPUT_DATA_1 for the construction of the joint histogram. 
%       Values that fall outside of the range are histogrammed in the first 
%       and last bins. If specified as input parameters, the same value is
%       return as output. 
%
%    HIST_LIMITS_2: [min_2,max_2] minimum and maximum values of variable 
%       INPUT_DATA_2 for the construction of the joint histogram. 
%       Values that fall outside of the range are histogrammed in the first 
%       and last bins. If specified as input parameters, the same value is
%       return as output. 
%
%    GRADIENT_JOINT_ENTROPY_2: gradient of the Joint Entropy with respect 
%       to INPUT_DATA_2. GRADIENT_JOINT_ENTROPY_2 has the same size as
%       INPUT_DATA_2 (and INPUT_DATA_1). 
%
%    ENTROPY_1: (marginal) Entropy of INPUT_DATA_1 calculated on the (marginal) kernel
%    density estimate. 
%
%    ENTROPY_2: (marginal) Entropy of INPUT_DATA_2 calculated on the (marginal) kernel
%    density estimate. 
%
%    GRADIENT_ENTROPY_1: gradient of the (marginal) Entropy of INPUT_DATA_1 with respect 
%       to INPUT_DATA_1. GRADIENT_ENTROPY_1 has the same size as INPUT_DATA_1 (and INPUT_DATA_2). 
%
%    GRADIENT_ENTROPY_2: gradient of the (marginal) Entropy of INPUT_DATA_2 with respect 
%       to INPUT_DATA_2. GRADIENT_ENTROPY_1 has the same size as
%       INPUT_DATA_2 (and INPUT_DATA_1). 
%
%
%
%Algorithm notes
%    Efficient estimation of the bivariate density and of the gradient of 
%    the Joint Entropy based on [1] and [2]. 
%    Automatic bandwidth selection based on [1]. 
%
%Reference
%    [1] Z. I. Botev, J. F. Grotowski and D. P. Kroese, 
%    "Kernel Density Estimation via Diffusion", Annals of Statistics, 2009.
%    
%    [2] S. Shwartz, M. Zibulevsky, Y.Y. Schechner,
%    "Fast kernel entropy estimation and optimization", Signal Processing,
%    5(85), 1045--1058, 2005. 
%
%    [3] S. Pedemonte, A. Bousse, K. Erlandsson, M. Modat, S. Arridge, B. Hutton, S. Ourselin, 
%    "GPU Accelerated Rotation-Based Emission Tomography Reconstruction", NSS/MIC 2010
%
%Example
%    N       = 128;
%    X1 = et_spherical_phantom(N,N,N,N/8,1,0.1,N/2,N/2,N/2)+rand(N,N,N);
%    X2 = et_spherical_phantom(N,N,N,N/4,2,0.2,N/2,N/2,N/2)+rand(N,N,N);
%    [JE,gradient_JE,density,histogram] = et_joint_entropy_nonparametric(X1,X2); 
%    figure; imagesc(gradient_JE(:,:,N/2)); colormap gray; 
%    figure; imagesc((density-min(density(:))).^0.4); colormap gray; 
%    figure; plot(sum(density,1),'r'); hold on; plot(sum(density,2),'r'); 
%    plot(sum(histogram,1)); hold on; plot(sum(histogram,2)); 
%
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


    default_N_bins = 256; 
    default_bandwidth = 0.1; 
    auto_bandwidth = 1; 
    auto_hist_limits = 1; 
    hist_limits_auto_range = 0.25; 
    
    if nargin >= 3 
        N_bins_1 = varargin{1}; 
    else
        N_bins_1 = default_N_bins; 
    end
    if nargin >= 4 
        N_bins_2 = varargin{2}; 
    else
        N_bins_2 = default_N_bins; 
    end
    if nargin >= 6 
        bandwidth_1 = varargin{3}; 
        bandwidth_2 = varargin{4}; 
        auto_bandwidth = 0; 
    else
        auto_bandwidth = 1; 
    end    
    if nargin >= 8
        hist_limits_1 = varargin{5}; 
        hist_limits_2 = varargin{6}; 
        auto_hist_limits = 0;        
    else
        auto_hist_limits = 1; 
    end

    % Make sure that input data and parameters are coherent
    if size(input_data_1) ~= size(input_data_2)
        error('Size of the two input arrays must be the same.');
    end
    input_shape = size(input_data_1); 
    input_data_1 = input_data_1(:);
    input_data_2 = input_data_2(:);
    N = length(input_data_1); 
    
    if N_bins_1 ~= 2^ceil(log2(N_bins_1))
        N_bins_1 = 2^ceil(log2(N_bins_1)); 
        warning(sprintf('Number of bins N_bins_1 must be a power of 2, using %d',N_bins_1));
    end
    if N_bins_2 ~= 2^ceil(log2(N_bins_2))
        N_bins_2 = 2^ceil(log2(N_bins_2)); 
        warning(sprintf('Number of bins N_bins_2 must be a power of 2, using %d',N_bins_2));
    end    
    
    % Calculate joint histogram
    limits_1 = [min(input_data_1),max(input_data_1)]; 
    limits_2 = [min(input_data_2),max(input_data_2)]; 
    range_1 = limits_1(2) - limits_1(1); 
    range_2 = limits_2(2) - limits_2(1); 
    if auto_hist_limits
       hist_limits_1 = [min(input_data_1)-hist_limits_auto_range*range_1,max(input_data_1)+hist_limits_auto_range*range_1]; 
       hist_limits_2 = [min(input_data_2)-hist_limits_auto_range*range_2,max(input_data_2)+hist_limits_auto_range*range_2]; 
    end
    
    discrete_data_1 = floor( (input_data_1-hist_limits_1(1)+10*eps)*(N_bins_1-10*eps)/(hist_limits_1(2)-hist_limits_1(1)) );
    discrete_data_2 = floor( (input_data_2-hist_limits_2(1)+10*eps)*(N_bins_2-10*eps)/(hist_limits_2(2)-hist_limits_2(1)) ); 
    discrete_data_1(discrete_data_1>=(N_bins_1-1)) = N_bins_1-1; 
    discrete_data_1 = discrete_data_1+1; 
    discrete_data_2(discrete_data_2>=(N_bins_2-1)) = N_bins_2-1; 
    discrete_data_2 = discrete_data_2+1; 

    discrete_data_array = discrete_data_1 + N_bins_1*(discrete_data_2-1); 
    joint_histogram_array = histc(discrete_data_array,[1:1:N_bins_1*N_bins_2]); 
    joint_histogram = reshape(joint_histogram_array,N_bins_1,N_bins_2); 

    % Automatic bandwidth 
    if auto_bandwidth 
        % discrete cosine transform of initial data 
        sigma_1 = N_bins_1 * default_bandwidth / 2.35482; 
        sigma_2 = N_bins_2 * default_bandwidth / 2.35482; 
        warning(sprintf('Automatic adaptive bandwidth not implemented, using default bandwidth: %f.\n',default_bandwidth)); 
    else
        sigma_1 = N_bins_1 * bandwidth_1 / 2.35482;
        sigma_2 = N_bins_2 * bandwidth_2 / 2.35482;
    end
    
    % Kernel density estimate
    size_kernel_1 = round(sigma_1*7.2); 
    if mod(size_kernel_1,2)==0
        size_kernel_1 = size_kernel_1+1;
    end
    size_kernel_2 = round(sigma_2*7.2); 
    if mod(size_kernel_2,2)==0
        size_kernel_2 = size_kernel_2+1;
    end
    kernel_support_1 = [-(size_kernel_1-1)/2:1:(size_kernel_1-1)/2]; 
    kernel_support_2 = [-(size_kernel_2-1)/2:1:(size_kernel_2-1)/2]; 
    kernel_1 = exp(-0.5*kernel_support_1.^2/(sigma_1^2)); 
    kernel_1 = kernel_1/sum(kernel_1(:)); 
    kernel_2 = exp(-0.5*kernel_support_2.^2/(sigma_2^2)); 
    kernel_2 = kernel_2/sum(kernel_2(:)); 

    joint_density = convnsep({kernel_1,kernel_2},joint_histogram,'same'); 
    joint_density = joint_density / sum(joint_density(:)); 
   
    % Estimate joint entropy 
    log_joint_density = zeros(size(joint_density));
    mask_log = joint_density>eps;
    log_joint_density(mask_log) = log(joint_density(mask_log)); 
    joint_entropy = -sum(joint_density(:) .* log_joint_density(:) ); 

    % Compute the gradient of the Joint Entropy with respect to data 1
    if nargout >= 2
        kernel_diff_1 = exp(-0.5*kernel_support_1.^2/(sigma_1^2)); 
        kernel_diff_1 = - kernel_support_1 .* kernel_diff_1/sum(kernel_diff_1(:)); 
        kernel_diff_2 = exp(-0.5*kernel_support_2.^2/(sigma_2^2)); 
        kernel_diff_2 = kernel_diff_2/sum(kernel_diff_2(:)); 
        gradient = convnsep({kernel_diff_1,kernel_diff_2}, 1+log_joint_density, 'same'); 
        gradient_joint_entropy_1 = reshape(gradient( discrete_data_array ),input_shape); 
    end    

    % Compute the gradient of the Joint Entropy with respect to data 2
    if nargout >= 7
        kernel_diff_1 = exp(-0.5*kernel_support_1.^2/(sigma_1^2)); 
        kernel_diff_1 = kernel_diff_1/sum(kernel_diff_1(:)); 
        kernel_diff_2 = exp(-0.5*kernel_support_2.^2/(sigma_2^2)); 
        kernel_diff_2 = - kernel_support_2 .* kernel_diff_2/sum(kernel_diff_2(:)); 
        gradient = convnsep({kernel_diff_1,kernel_diff_2}, 1+log_joint_density, 'same'); 
        gradient_joint_entropy_2 = reshape(gradient( discrete_data_array ),input_shape); 
    end

    % Compute the Entropy of data 1 and data 2
    if nargout >= 8
        density_1 = sum(joint_density,2); density_1 = density_1/sum(density_1(:));
        log_density_1 = zeros(size(density_1));
        mask_log = density_1>eps;
        log_density_1(mask_log) = log(density_1(mask_log)); 
        entropy_1 = -sum(density_1(:) .* log_density_1(:) ); 
    end
    if nargout >= 9    
        density_2 = sum(joint_density,1); density_2 = density_2/sum(density_2(:));
        log_density_2 = zeros(size(density_2));
        mask_log = density_2>eps;
        log_density_2(mask_log) = log(density_2(mask_log)); 
        entropy_2 = -sum(density_2(:) .* log_density_2(:) ); 
    end

    % Compute the gradient of the Entropy with respect to data 1 and data 2    
    if nargout >= 10
        kernel_diff_1 = exp(-0.5*kernel_support_1.^2/(sigma_1^2)); 
        kernel_diff_1 = - kernel_support_1 .* kernel_diff_1/sum(kernel_diff_1(:)); 
        gradient = conv(1+log_density_1, kernel_diff_1, 'same'); 
        gradient_entropy_1 = reshape(gradient( discrete_data_1 ),input_shape); 
    end
    if nargout >= 11    
        kernel_diff_2 = exp(-0.5*kernel_support_2.^2/(sigma_2^2)); 
        kernel_diff_2 = - kernel_support_2 .* kernel_diff_2/sum(kernel_diff_2(:)); 
        gradient = conv(1+log_density_2, kernel_diff_2, 'same'); 
        gradient_entropy_2 = reshape(gradient( discrete_data_2 ),input_shape);         
    end
    
    % Optional outputs
    if nargout >= 2
        varargout{1} = gradient_joint_entropy_1;
    end
    
    if nargout >= 3
        varargout{2} = joint_density;
    end
    
    if nargout >= 4
        varargout{3} = joint_histogram; 
    end
    
    if nargout >= 5
        varargout{4} = hist_limits_1; 
    end
    
    if nargout >= 6
        varargout{5} = hist_limits_2; 
    end
    
    if nargout >= 7
        varargout{6} = gradient_joint_entropy_2; 
    end

    if nargout >= 8
        varargout{7} = entropy_1; 
    end
    
    if nargout >= 9
        varargout{8} = entropy_2; 
    end    
    
    if nargout >= 10
        varargout{9} = gradient_entropy_1;
    end

    if nargout >= 11
        varargout{10} = gradient_entropy_2;
    end    
end






%#######################################
function J = convnsep(h,V,type)
lh=length(h);

%input validation
assert(lh==ndims(V),'The number of kernels does not match the array dimensionality.');

L=nan(1,lh);
for j=1:lh,
    L(j)=(length(h{j})-1)/2;
end
V=padarray(V,L);
J = convnsepsame(h,V);

%implicit behaviour: if no 'type' input, then type=='full' (don't trim the
%result)
if nargin>2
    switch type
        case 'full'
            %do nothing
        case 'valid'
            J=trimarray(J,L*2);
        case 'same'
            J=trimarray(J,L);
        otherwise
    end
end

end

%Perform convolution while keeping the array size the same (i.e. discarding
%boundary samples)
function J = convnsepsame(h,V)
J=V;
sz=size(V);
n=length(sz);
indx2=nan(1,n);
for k=1:n
    %dimensions other k-th dimension, along which convolution will happen:
    otherdims = 1:n; otherdims(k)=[];
    % permute order: place k-th dimension as 1st, followed by all others:
    indx1=[k otherdims];
    % inverse permute order:
    indx2(indx1)=1:n;
    %perform convolution along k-th dimension:
    %
    %1. permute dimensions to place k-th dimension as 1st
    J = permute(J,indx1);
    %2. create a 2D array (i.e. "stack" all other dimensions, other than
    %k-th:
    J = reshape(J,sz(k),prod(sz(otherdims)));
    %3. perform 2D convolution with k-th kernel along the first dimension
    J = conv2(h{k},1,J,'same');
    %4. undo the "flattening" of step 2
    J = reshape(J,sz(indx1));
    %5. undo the permutation of step 1.
    J = permute(J,indx2);
end
end

%extract only the central portion of the array V, based on kernels whose
%lengths are defined in L
function V = trimarray(V,L)
str='';
for j=1:ndims(V)
    str=[str num2str(L(j)+1) ':' num2str(size(V,j)-L(j)) ','];
end
str=str(1:end-1); %remove last coma

eval(['V=V(' str ');']);
end


