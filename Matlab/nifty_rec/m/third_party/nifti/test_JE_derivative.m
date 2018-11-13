function test_JE_derivative

clc;

image1 = 1 + 63*rand(10,10);
image2 = 1 + 63*rand(10,10);


figure;imagesc(image1);colormap gray; title('Reference');
figure;imagesc(image2);colormap gray; title('Initial reconstructed');

for i=1:20;
    [JE_scalar, JE_gradient] = JE_derivative(image1,image2);
    fprintf('JE value = %g\n', JE_scalar);
    image2= image2 - JE_gradient;
end

figure;imagesc(image2);colormap gray; title('Final reconstructed');