
function density = et_normal_density(mu, sigma, N_points, min_value, max_value)

x = linspace(min_value,max_value,N_points); 
density = (1/sqrt(2*pi*sigma^2))*exp(-0.5*(x-mu).^2/(sigma^2)) *(max_value-min_value)/N_points; 

