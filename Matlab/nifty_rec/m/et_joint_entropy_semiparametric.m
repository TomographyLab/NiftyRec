function [joint_entropy, varargout] = et_joint_entropy_semiparametric(input_data_NONPARAM, input_data_PARAM, prior_labels, mu, sigma, varargin) 
%ET_JOINT_ENTROPY_SEMIPARAMETRIC 
%    Semiparametric kernel density estimator with diagonal covariance
%    matrix. 
%    The function returns the Joint Entropy between the two variables measured 
%    on the semiparametric density estimate; the gradient of the Joint Entropy with
%    respect of the input variable INPUT_DATA_PARAM and, optionally, 
%    the bivariate density estimate. 
%    The function also optionally returns the marginal densities, the 
%    marginal entropies and the gradient of the marginal entropies with
%    respect to each of the input variables. 
%
%Description
%
%    [JOINT_ENTROPY, GRADIENT_JOINT_ENTROPY_NONPARAM, JOINT_DENSITY, JOINT_HISTOGRAM, 
%    HIST_LIMITS_NONPARAM, HIST_LIMITS_PARAM, GRADIENT_JOINT_ENTROPY_PARAM, ENTROPY_NONPARAM, ENTROPY_PARAM, 
%    GRADIENT_ENTROPY_NONPARAM, GRADIENT_ENTROPY_PARAM] = et_jointentropy_semiparametric(INPUT_DATA_NONPARAM, 
%    INPUT_DATA_PARAM, LABELS_PROBABILITY, N_BINS_NONPARAM, N_BINS_PARAM, BANDWIDTH, HIST_LIMITS_NONPARAM, HIST_LIMITS_PARAM); 
%
%  All input arguments are optional except for the input data INPUT_DATA_NONPARAM
%  and INPUT_DATA_PARAM. 
%  All output arguments are optional except for the first (JOINT ENTROPY). 
%
%
%  Input arguments: 
%
%    INPUT_DATA_NONPARAM: this is a [Nx1] array of continuous values. 
%       INPUT_DATA_NONPARAM can also be specified as a multi-dimensional array (e.g. 2D image,
%       3D image). In this case the gradient (output) will
%       have the same size as the input data. The size of INPUT_DATA_NONPARAM must 
%       be equal to the size of INPUT_DATA_PARAM. 
%
%    INPUT_DATA_PARAM: this should be a [Nx1] array of continuous values.
%       INPUT_DATA_PARAM can also be specified as a multi-dimensional array (e.g. 2D image,
%       3D image). In this case the gradient, returned by the function, will
%       have the same size as the input data. The size of INPUT_DATA_PARAM must 
%       be equal to the sizeof INPUT_DATA_NONPARAM. 
%
%    LABELS_PROBABILITY: probability that hidden label is in state 1-of-N_classes.
%
%
%  Optional input arguments: 
%
%    N_BINS_NONPARAM, N_BINS_PARAM: number of bins for the discretisation of the joint
%       histogram, along the axis INPUT_DATA_NONPARAM and INPUT_DATA_PARAM. 
%       The joint histogram is constructed on a square grid of
%       N_BINS_NONPARAM*N_BINS_PARAM elements. 
%       Construction of the joint histogram and smoothing are a computationally
%       efficient way to discretise the kernel density estimation problem. 
%       The higher the value of N_BINS_NONPARAM and N_BINS_PARAM, the more accurate the estimate of the
%       joint density. The lower the value of N_BINS_NONPARAM and N_BINS_PARAM, the faster the
%       algorithm and the lower the memory requirements. 
%       Default value: 256. Typical values: 256 or 512. 
%
%    BANDWIDTH: bandwidth (Full Width at Half Maximum - FWHM) of the 
%       Gaussian kernel for variable INPUT_DATA_PARAM. 
%       FWHM=1 means that the FWHM of the kernel is as large as the domain of
%       the histogram. 
%       Typical value: 0.02 - 0.2 
%       If not set, then the bandwidth is selected automatically. 
%
%    HIST_LIMITS_NONPARAM: [min_NONPARAM,max_NONPARAM] minimum and maximum values of variable 
%       INPUT_DATA_NONPARAM for the construction of the joint histogram. 
%       They are autimatically set if not specified. Values that fall outside 
%       of the range are histogrammed in the first and last bins. 
%
%    HIST_LIMITS_PARAM: [min_PARAM,max_PARAM] minimum and maximum values of variable 
%       INPUT_DATA_PARAM for the construction of the joint histogram. 
%       They are autimatically set if not specified. Values that fall outside 
%       of the range are histogrammed in the first and last bins.  
%
%
%  Output arguments: 
%
%    JOINT_ENTROPY Joint Entropy of INPUT_DATA_NONPARAM and INPUT_DATA_PARAM, 
%       calculated according to the kernel joint density estimate. 
%       It is a scalar value. 
% 
%
%  Optional output arguments: 
%
%    GRADIENT_JOINT_ENTROPY_NONPARAM gradient of the Joint Entropy with respect 
%       to INPUT_DATA_NONPARAM. GRADIENT_JOINT_ENTROPY_NONPARAM has the same size as
%       INPUT_DATA_NONPARAM (and INPUT_DATA_PARAM). 
%
%    JOINT_DENSITY: kernel density estimate of (INPUT_DATA_NONPARAM,INPUT_DATA_PARAM).
%
%    JOINT_HISTOGRAM: joint histogram of (INPUT_DATA_NONPARAM,INPUT_DATA_PARAM). 
% 
%    HIST_LIMITS_NONPARAM: [min_NONPARAM,max_NONPARAM] minimum and maximum values of variable 
%       INPUT_DATA_NONPARAM for the construction of the joint histogram. 
%       Values that fall outside of the range are histogrammed in the first 
%       and last bins. If specified as input parameters, the same value is
%       return as output. 
%
%    HIST_LIMITS_PARAM: [min_PARAM,max_PARAM] minimum and maximum values of variable 
%       INPUT_DATA_PARAM for the construction of the joint histogram. 
%       Values that fall outside of the range are histogrammed in the first 
%       and last bins. If specified as input parameters, the same value is
%       return as output. 
%
%    GRADIENT_JOINT_ENTROPY_PARAM: gradient of the Joint Entropy with respect 
%       to INPUT_DATA_PARAM. GRADIENT_JOINT_ENTROPY_PARAM has the same size as
%       INPUT_DATA_PARAM (and INPUT_DATA_NONPARAM). 
%
%    ENTROPY_NONPARAM: (marginal) Entropy of INPUT_DATA_NONPARAM calculated on the (marginal) kernel
%    density estimate. 
%
%    ENTROPY_PARAM: (marginal) Entropy of INPUT_DATA_PARAM calculated on the (marginal) kernel
%    density estimate. 
%
%    GRADIENT_ENTROPY_NONPARAM: gradient of the (marginal) Entropy of INPUT_DATA_NONPARAM with respect 
%       to INPUT_DATA_NONPARAM. GRADIENT_ENTROPY_NONPARAM has the same size as INPUT_DATA_NONPARAM (and INPUT_DATA_PARAM). 
%
%    GRADIENT_ENTROPY_PARAM: gradient of the (marginal) Entropy of INPUT_DATA_PARAM with respect 
%       to INPUT_DATA_PARAM. GRADIENT_ENTROPY_NONPARAM has the same size as
%       INPUT_DATA_PARAM (and INPUT_DATA_NONPARAM). 
%
%
%
%Algorithm notes
%
%    Efficient estimation of the bivariate density and of the gradient of 
%    the Joint Entropy based on [1] and [2]. 
%    Automatic bandwidth selection based on [1]. 
%
%Reference
%
%    [3] S. Pedemonte, A. Bousse, K. Erlandsson, M. Modat, S. Arridge, B. Hutton, S. Ourselin, 
%    "GPU Accelerated Rotation-Based Emission Tomography Reconstruction", NSS/MIC 2010
%
%
%
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK

    default_N_bins = 256; 
    default_bandwidth = 0.1; 
    auto_bandwidth = 1; 
    auto_hist_limits = 1; 
    hist_limits_auto_range = 0.25; 
    
    if nargin >= 5
        N_bins_NONPARAM = varargin{1}; 
    else
        N_bins_NONPARAM = default_N_bins; 
    end
    if nargin >= 6 
        N_bins_PARAM = varargin{2}; 
    else
        N_bins_PARAM = default_N_bins; 
    end
    if nargin >= 8 
        bandwidth_NONPARAM = varargin{3}; 
        auto_bandwidth = 0; 
    else
        auto_bandwidth = 1; 
    end    
    if nargin >= 9
        hist_limits_NONPARAM = varargin{4}; 
        hist_limits_PARAM = varargin{5}; 
        auto_hist_limits = 0;        
    else
        auto_hist_limits = 1; 
    end
    if nargin >= 11
        class_mask = varargin{6}; 
    else
        class_mask = 0;
    end
        
    % Make sure that input data and parameters are coherent
    if size(input_data_NONPARAM) ~= size(input_data_PARAM)
        error('Size of the two input arrays must be the same.');
    end
    input_shape = size(input_data_NONPARAM); 
    input_data_NONPARAM = input_data_NONPARAM(:);
    input_data_PARAM = input_data_PARAM(:);
    
    if N_bins_NONPARAM ~= 2^ceil(log2(N_bins_NONPARAM))
        N_bins_NONPARAM = 2^ceil(log2(N_bins_NONPARAM)); 
        warning(sprintf('Number of bins N_bins_NONPARAM must be a power of 2, using %d',N_bins_NONPARAM));
    end
    if N_bins_PARAM ~= 2^ceil(log2(N_bins_PARAM))
        N_bins_PARAM = 2^ceil(log2(N_bins_PARAM)); 
        warning(sprintf('Number of bins N_bins_PARAM must be a power of 2, using %d',N_bins_PARAM));
    end    
    N_classes = size(prior_labels,length(size(prior_labels))); 

    % Automatic bandwidth 
    if auto_bandwidth 
        % discrete cosine transform of initial data 
        sigma_NONPARAM = N_bins_NONPARAM * default_bandwidth / 2.35482; 
        warning(sprintf('Automatic adaptive bandwidth not implemented, using default bandwidth: %f.\n',default_bandwidth)); 
    else
        sigma_NONPARAM = N_bins_NONPARAM * bandwidth_NONPARAM / 2.35482;
    end
    
    % Calculate class membership 
    %membership = zeros([size(input_data_PARAM),N_classes]); 
    membership = prior_labels;  % !!!!
    
    membership = reshape(membership,prod(input_shape(:)),N_classes);
    mixing_coefficient = zeros(1,N_classes); 
    for class=1:N_classes
        mixing_coefficient(class)=sum(membership(:,class))/sum(membership(:));
    end
    
    % Estimate the conditional probability distribution nonparametric   p(NONPARAM|class)
    if auto_hist_limits
       hist_limits_NONPARAM = [min(input_data_NONPARAM)-hist_limits_auto_range*range_NONPARAM,max(input_data_NONPARAM)+hist_limits_auto_range*range_NONPARAM]; 
       hist_limits_PARAM = [min(input_data_PARAM)-hist_limits_auto_range*range_PARAM,max(input_data_PARAM)+hist_limits_auto_range*range_PARAM]; 
    end    
    
    histograms_NONPARAM = et_histogram_weighted_mex(input_data_NONPARAM(:),membership,N_bins_NONPARAM,hist_limits_NONPARAM(1),hist_limits_NONPARAM(2)); 

    size_kernel_NONPARAM = round(sigma_NONPARAM*7.2); 
    if mod(size_kernel_NONPARAM,2)==0
        size_kernel_NONPARAM = size_kernel_NONPARAM+1;
    end
    kernel_support_NONPARAM = [-(size_kernel_NONPARAM-1)/2:1:(size_kernel_NONPARAM-1)/2]; 
    kernel_NONPARAM = exp(-0.5*kernel_support_NONPARAM.^2/(sigma_NONPARAM^2)); 
    kernel_NONPARAM = kernel_NONPARAM/sum(kernel_NONPARAM(:));     
    
    density_NONPARAM = zeros([N_bins_NONPARAM,N_classes]);
    for class=1:N_classes 
        density_NONPARAM(:,class) = conv(histograms_NONPARAM(:,class),kernel_NONPARAM,'same'); 
        density_NONPARAM(:,class) = density_NONPARAM(:,class) / sum(density_NONPARAM(:,class));
    end

    epsilon=1e-8;
    discrete_data_NONPARAM = floor( (input_data_NONPARAM-hist_limits_NONPARAM(1)+10*epsilon)*(N_bins_NONPARAM-10*epsilon)/(hist_limits_NONPARAM(2)-hist_limits_NONPARAM(1)) );
    discrete_data_PARAM = floor( (input_data_PARAM-hist_limits_PARAM(1)+10*epsilon)*(N_bins_PARAM-10*epsilon)/(hist_limits_PARAM(2)-hist_limits_PARAM(1)) ); 
    discrete_data_NONPARAM(discrete_data_NONPARAM>=(N_bins_NONPARAM-1)) = N_bins_NONPARAM-1; 
    discrete_data_NONPARAM = discrete_data_NONPARAM+1; 
    discrete_data_PARAM(discrete_data_PARAM>=(N_bins_PARAM-1)) = N_bins_PARAM-1; 
    discrete_data_PARAM = discrete_data_PARAM+1; 
    discrete_data_array = discrete_data_NONPARAM + N_bins_NONPARAM*(discrete_data_PARAM-1); 

    % Estimate the conditional probability distribution parametric   p(PARAM|class)
    density_PARAM = zeros(N_bins_PARAM,N_classes); 
    for class=1:N_classes
        density_PARAM(:,class) = et_normal_density(mu(class),sigma(class),N_bins_PARAM,hist_limits_PARAM(1),hist_limits_PARAM(2)); 
        density_PARAM(:,class) = density_PARAM(:,class) / sum(density_PARAM(:,class));
    end
        
    % Estimate the joint conditional probability distribution   p(NONPARAM,PARAM|class)
    density_joint_conditional = zeros(N_bins_NONPARAM,N_bins_PARAM,N_classes); 
    for class=1:N_classes
        density_joint_conditional(:,:,class) = density_NONPARAM(:,class)*density_PARAM(:,class)'; 
        density_joint_conditional(:,:,class) = density_joint_conditional(:,:,class) / sum(sum(density_joint_conditional(:,:,class)));
    end
        
    % Estimate the joint probability distribution   p(NONPARAM,PARAM)
    density_joint = zeros(N_bins_NONPARAM,N_bins_PARAM); 
    for class=1:N_classes
        density_joint(:,:) = density_joint(:,:) + mixing_coefficient(class)*density_joint_conditional(:,:,class); 
    end
    density_joint = density_joint / sum(density_joint(:));
    
    % Estimate the joint entropy 
    log_density_joint = zeros(size(density_joint)) + log(eps);
    mask_log = density_joint>eps;
    log_density_joint(mask_log) = log(density_joint(mask_log)); 
    joint_entropy = -sum(density_joint(:) .* log_density_joint(:) ); 
    
    % Compute the gradient of the Joint Entropy with respect to input_data_NONPARAM
    if nargout >= 2
        % differential kernel 
        kernel_diff_NONPARAM = exp(-0.5*kernel_support_NONPARAM.^2/(sigma_NONPARAM^2)); 
        kernel_diff_NONPARAM = - kernel_support_NONPARAM .* kernel_diff_NONPARAM/sum(kernel_diff_NONPARAM(:));                     
        % compute k-th additive term of the gradient of the joint entropy
        gradients_class = zeros(N_bins_NONPARAM,N_bins_PARAM,N_classes); 
        for class=1:N_classes
            % multiply by conditional density of PARAM variable
            gradients_class(:,:,class) = (1+log_density_joint) ;%.* repmat(density_PARAM(:,class)',N_bins_NONPARAM,1); 
            % convolve with differential kernel
            gradients_class(:,:,class) = convnsep({kernel_diff_NONPARAM,[0,1,0]},reshape(gradients_class(:,:,class),N_bins_NONPARAM,N_bins_PARAM),'same'); 
        end
        % back-propagate the gradient
        gradients_joint_entropy_NONPARAM_class = zeros([prod(input_shape),N_classes]); 
        for class=1:N_classes
            gradient = gradients_class(:,:,class); 
            gradients_joint_entropy_NONPARAM_class(:,class) = gradient( discrete_data_array ); 
        end

        % combine the N_classes terms of the gradient of the joint entrop (combine them in the space of the input data) 
        gradient_joint_entropy_NONPARAM = zeros(size(input_data_NONPARAM)); 
        if class_mask==0
            classes = 1:N_classes;
        else
            classes = class_mask;
        end
        for class=classes
                gradient_joint_entropy_NONPARAM = gradient_joint_entropy_NONPARAM + mixing_coefficient(class) .* membership(:,class) .* gradients_joint_entropy_NONPARAM_class(:,class); 
        end
        % reshape as input
        membership = reshape(membership,[input_shape,N_classes]); 
        gradients_joint_entropy_NONPARAM_class = reshape(gradients_joint_entropy_NONPARAM_class,[input_shape,N_classes]); 
        gradient_joint_entropy_NONPARAM = reshape(gradient_joint_entropy_NONPARAM,input_shape); 
    end 
    
 
    % Optional outputs
    if nargout >= 2
        varargout{1} = gradient_joint_entropy_NONPARAM;
    end
    
    if nargout >= 3
        varargout{2} = density_joint;
    end
    
    if nargout >= 4
        varargout{3} = histograms_NONPARAM; 
    end
    
    if nargout >= 5
        varargout{4} = hist_limits_NONPARAM; 
    end
    
    if nargout >= 6
        varargout{5} = hist_limits_PARAM; 
    end

    if nargout >= 7
        varargout{6} = gradients_joint_entropy_NONPARAM_class; 
    end

    if nargout >= 8
        varargout{7} = density_joint_conditional; 
    end    
    
    if nargout >= 9
        varargout{10} = membership; 
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

