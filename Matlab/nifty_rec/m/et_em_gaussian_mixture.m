function [mu,sigma,mix,membership] = et_em_gaussian_mixture(image, n_classes, em_steps, mu, sigma)

data = image(:);

data_size = numel(data);
if nargin < 4 
    m = mean(data)-3*std(data);
    M = mean(data)+3*std(data);
    delta = (M-m)/n_classes;
    mu = (m+delta/2:delta:M-delta/4);
    sigma = ones(1,n_classes) * std(data(:)) / n_classes;
end
mix = ones(1,n_classes)/n_classes;
membership = zeros(n_classes,data_size);
pdf = zeros(1,data_size); 

for step = 1:em_steps
    fprintf('\nStep %d',step);
    %E step
    pdf = pdf * 0;
    for class = 1:n_classes
        membership(class,:) = mix(class) * (normpdf(data,mu(class),sigma(class)) + eps);
        pdf = pdf + membership(class,:);
    end
    for class = 1:n_classes
        membership(class,:) = membership(class,:) ./ pdf;
    end
    %M step
    for class = 1:n_classes
        mix(class) = (1/data_size) .* sum(membership(class,:));
        sigma(class)  = ((1/data_size) .* (1./mix(class)) .* sum(reshape(membership(class,:),data_size,1).*(data-mu(class)).^2)).^0.5;
        mu(class)  = (1/data_size) .* (1./mix(class)) .* sum(reshape(membership(class,:),data_size,1).*data);
    end
end

membership_array = membership;
membership = cell(1,n_classes);
for i = 1:n_classes
    membership{i} = reshape(membership_array(i,:),size(image));
end

fprintf('\nDone\n');
end