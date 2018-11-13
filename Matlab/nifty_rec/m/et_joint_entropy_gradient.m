function [joint_entropy, gradient_joint_entropy_X2, varargout] = et_joint_entropy_gradient(X1, X2, N_bins, bandwidth_X1, bandwidth_X2, enable_automatic_bandwidth, min_grid_histogram, max_grid_histogram)
%ET_JOINT_ENTROPY_GRADIENT
%    Fast and accurate bivariate kernel density estimator with Gaussian kernel and diagonal
%    bandwidth matrix. The function returns the Joint Entropy between the two variables measured 
%    on the kernel density estimate; the gradient of the Joint Entropy with
%    respect of variable 'floating' and, optionally, the bivariate density estimate. 
%    The two bandwidth parameters are optionally estimated automatically with the
%    algorithm described in Botev 2009 [1]. 
%
%Description
%    [JOINT_ENTROPY, GRADIENT_JOINT_ENTROPY_X2, DENSITY, BANDWIDTH_X1, BANDWIDTH_X2, 
%    HISTOGRAM, GRID] = et_joint_entropy_gradient(X1, X2, N_BINS, BANDWIDTH_X1, 
%    BANDWIDTH_X2, ENABLE_AUTOMATIC_BANDWIDTH)
%
%
%  Input arguments: 
%
%    X1 input data 1. X1 should be a [Nx1] array of continuous values.
%       X1 can also be specified as a multi-dimensional array (e.g. 2D image,
%       3D image). In this case the gradient, returned by the function, will
%       have the same size as the input data. The size of X1 must 
%       be equal to the size of X2. 
%
%    X2 input data 2. X2 should be a [Nx1] array of continuous values.
%       X2 can also be specified as a multi-dimensional array (e.g. 2D image,
%       3D image). In this case the gradient, returned by the function, will
%       have the same size as the input data. The size of X2 must 
%       be equal to the sizeof X1. 
%       The function returns the derivatives of the Joint Entropy with respect
%       to the elements of X2. 
%
%    N_BINS number of bins for the discretisation of the joint histogram. 
%       The joint histogram is constructed on a square grid of N_BINS*N_BINS
%       elements. 
%       Construction of the joint histogram and smoothing are a computationally
%       efficient way to discretise the kernel density estimation problem. 
%       The higher the value of N_BINS, the more accurate the estimates of the
%       joint density. The lower the value of N_BINS, the faster the
%       algorithm. 
%       Typical values: 256, 512. 
%
%    BANDWIDTH_X1 kernel bandwidth for variable X1. Not used if
%       ENABLE_OPTIMAL_BANDWIDTH not equal 0;
%
%    BANDWIDTH_X2 kernel bandwidth for variable X2. Not used if
%       ENABLE_OPTIMAL_BANDWIDTH not equal 0;
%
%    ENABLE_AUTOMATIC_BANDWIDTH flag to enable/disable automatic bandwidth 
%       estimation. If nonzero, the bandwidths of the kernels are estimated
%       automatically with the algorithm described in [1], ignoring the
%       values of BANDWIDTH_X1 and BANDWIDTH_X2. 
%
%    MIN_GRID_HISTOGRAM [min_X1,min_X2] minimum values for the construction
%    of the joint histogram. They are autimatically set if not specified. 
%
%    MAX_GRID_HISTOGRAM [max_X1,max_X2] maximum values for the construction
%    of the joint histogram. They are autimatically set if not specified.  
%
%
%  Output arguments: 
%
%    JOINT_ENTROPY Joint Entropy of X1 and X2, calculated according to the
%       kernel joint density estimate. This is a scalar value. 
% 
%    GRADIENT_JOINT_ENTROPY_X2 gradient of the Joint Entropy with respect 
%       to X2. GRADIENT_JOINT_ENTROPY_X2 has the same size as X1
%       and X2. 
%
%
%  Optional output arguments: 
%
%    DENSITY kernel density estimate of (X1,X2)
%
%    HISTOGRAM joint histogram of (X1,X2), calculated on the 2D grid GRID. 
%
%    GRID_HISTOGRAM_X1 grid for the calculation of the joint histogram. 
%
%    GRID_HISTOGRAM_X2 grid for the calculation of the joint histogram. 
%
%    BANDWIDTH_X1 bandwidth of the kernel for X1, equal to the input value
%    if ENABLE_AUTOMATIC_BANDWIDTH==0. 
%
%    BANDWIDTH_X2 bandwidth of the kernel for X2, equal to the input value
%    if ENABLE_AUTOMATIC_BANDWIDTH==0. 
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
%    N_bins  = 100; 
%    X1 = et_spherical_phantom(N,N,N,N/8,1,0.1,N/2,N/2,N/2)+rand(N,N,N);
%    X2 = et_spherical_phantom(N,N,N,N/4,2,0.2,N/2,N/2,N/2)+rand(N,N,N);
%    [JE,gradient_JE,density,histogram] = et_joint_entropy_gradient(X1,X2,N_bins,10,10,0); 
%    figure; imagesc(gradient_JE(:,:,N/2)); colormap gray; 
%    figure; imagesc((density-min(density(:))).^0.4); colormap gray; 
%    figure; plot(sum(density,1),'r'); hold on; plot(sum(density,2),'r'); 
%    plot(sum(histogram,1)); hold on; plot(sum(histogram,2)); 
%
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


% verify consistency of the input
if not(size(X2)==size(X1))
    error('X1 and X2 must be of equal size.'); 
end

% round up n to the next power of 2;
N_bins = 2^ceil(log2(N_bins)); 

X = [X1(:),X2(:)]; 
N = size(X,1); 

% select automatically the grid for the joint histogram if not specified
if nargin<7
    MAX=max(X,[],1); MIN=min(X,[],1); MAX = MAX + 0.1; Range=MAX-MIN; 
    max_grid_histogram=MAX+Range/4; min_grid_histogram=MIN-Range/4;
end
scaling = max_grid_histogram-min_grid_histogram; 

%bin the data uniformly using regular grid;
transformed_X = (X-repmat(min_grid_histogram,N,1))./repmat(scaling,N,1); 
transformed_X_int = ceil( (X - repmat(min_grid_histogram,N,1))*N_bins./repmat(scaling,N,1) ); 
initial_X = ndhist(transformed_X,N_bins);

%optionally estimate kernel bandwidth automatically
if enable_automatic_bandwidth
    % discrete cosine transform of initial data 
    t_star=fzero( @(t)(t-evolve(t)),[0,0.1]);
    p_02=func([0,2],t_star);p_20=func([2,0],t_star); p_11=func([1,1],t_star);
    t_y=(p_02^(3/4)/(4*pi*N*p_20^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
    t_x=(p_20^(3/4)/(4*pi*N*p_02^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
    % smooth the discrete cosine transform of initial data using t_star
    %a_t=exp(-(0:n-1)'.^2*pi^2*t_x/2)*exp(-(0:n-1).^2*pi^2*t_y/2).*a;
    bandwidth_X1 = N_bins * sqrt(t_x);
    bandwidth_X2 = N_bins * sqrt(t_y);
end

%scaling
xx = (-N_bins/2:N_bins/2)';
yy = (-N_bins/2:N_bins/2);
%k = (1/sqrt(2*pi*bandwidth_X1^2)*exp(-xx.^2/(2*bandwidth_X1^2))) * (1/sqrt(2*pi*bandwidth_X2^2)*exp(-yy.^2/(2*bandwidth_X2^2)));
k = exp(-xx.^2/(2*bandwidth_X1^2)) * exp(-yy.^2/(2*bandwidth_X2^2));
k = k / sum(k(:));

density_fft = conv2fft(initial_X,k,'same');
density_fft = density_fft / sum(density_fft(:));

% apply the inverse discrete cosine transform
histogram = initial_X;
[X,Y]=meshgrid(min_grid_histogram(1):scaling(1)/(N_bins-1):max_grid_histogram(1),min_grid_histogram(2):scaling(2)/(N_bins-1):max_grid_histogram(2)); 

%compute gradient
%kd = - ( xx .* ((1/sqrt(2*pi*bandwidth_X1^2)*exp(-xx.^2/(2*bandwidth_X1^2)))) )  * (1/sqrt(2*pi*bandwidth_X2^2)*exp(-yy.^2/(2*bandwidth_X2^2)));
kd = - ( xx .* exp(-xx.^2/(2*bandwidth_X1^2)) ) * exp(-yy.^2/(2*bandwidth_X2^2));
p_tilde = zeros(size(density_fft));
mask = density_fft>0;
p_tilde(mask) = 1 + log(density_fft(mask));
gr = conv2fft(p_tilde,kd,'same');
mask = (transformed_X_int <= 0);
transformed_X_int(mask) = 1;
gradient_joint_entropy_X2 = gr( transformed_X_int(:,1) + N_bins*(transformed_X_int(:,2) ));
gradient_joint_entropy_X2 = reshape(gradient_joint_entropy_X2,size(X2)); 

%compute joint entropy
log_density=zeros(size(density_fft)); 
mask=density_fft>0; 
log_density(mask)=log(density_fft(mask)); 
joint_entropy = -sum(density_fft(:) .* log_density(:) ); 

%return optional variables 
varargout{1} = density_fft; 
varargout{2} = histogram; 
varargout{3} = [X]; 
varargout{4} = [Y]; 
varargout{5} = bandwidth_X1;
varargout{6} = bandwidth_X2;

end

%#######################################
function  [out,time]=evolve(t)
global N
Sum_func = func([0,2],t) + func([2,0],t) + 2*func([1,1],t);
time=(2*pi*N*Sum_func)^(-1/3);
out=(t-time)/time;
end

%#######################################
function out=func(s,t)
global N
if sum(s)<=4;
    Sum_func=func([s(1)+1,s(2)],t)+func([s(1),s(2)+1],t); const=(1+1/2^(sum(s)+1))/3;
    time=(-2*const*K(s(1))*K(s(2))/N/Sum_func)^(1/(2+sum(s)));
    out=psi(s,time);
else
    out=psi(s,t);
end

end

%#######################################
function out=psi(s,Time)
global I A2
% s is a vector
w=exp(-I*pi^2*Time).*[1,.5*ones(1,length(I)-1)];
wx=w.*(I.^s(1));
wy=w.*(I.^s(2));
out=(-1)^sum(s)*(wy*A2*wx')*pi^(2*sum(s));
end

%#######################################
function out=K(s)
out=(-1)^s*prod((1:2:2*s-1))/sqrt(2*pi);
end

%#######################################
function X=dct2d(X)
% computes the 2 dimensional discrete cosine transform of data
% data is an nd cube
[nrows,ncols]= size(X);
if nrows~=ncols
    error('data is not a square array!')
end
% Compute weights to multiply DFT coefficients
w = [1;2*(exp(-i*(1:nrows-1)*pi/(2*nrows))).'];
weight=w(:,ones(1,ncols));
X=dct1d(dct1d(X)')';
    function transform1d=dct1d(x)

        % Re-order the elements of the columns of x
        x = [ x(1:2:end,:); x(end:-2:2,:) ];

        % Multiply FFT by weights:
        transform1d = real(weight.* fft(x));
    end
end

%#######################################
function X = idct2d(X)
% computes the 2 dimensional inverse discrete cosine transform
[nrows,ncols]=size(X);
% Compute wieghts
w = exp(i*(0:nrows-1)*pi/(2*nrows)).';
weights=w(:,ones(1,ncols));
X=idct1d(idct1d(X)');
    function out=idct1d(x)
        y = real(ifft(weights.*x));
        out = zeros(nrows,ncols);
        out(1:2:nrows,:) = y(1:nrows/2,:);
        out(2:2:nrows,:) = y(nrows:-1:nrows/2+1,:);
    end
end

%#######################################
function binned_data=ndhist(data,M)
% this function computes the histogram
% of an n-dimensional data set;
% 'data' is nrows by n columns
% M is the number of bins used in each dimension
% so that 'binned_data' is a hypercube with
% size length equal to M;
[nrows,ncols]=size(data);
bins=zeros(nrows,ncols);
for i=1:ncols
    [dum,bins(:,i)] = histc(data(:,i),[0:1/M:1],1);
    bins(:,i) = min(bins(:,i),M);
end
% Combine the  vectors of 1D bin counts into a grid of nD bin
% counts.
binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));
end




















