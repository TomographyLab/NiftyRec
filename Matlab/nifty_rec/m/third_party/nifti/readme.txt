The zip file include:
- readme.txt  ---> this file
- conv2fft.m  ---> M-file for 2D FFT-based convolution
- convfft.m   ---> M-file for 1D FFT-based convolution


The convolution is performed in frequency domain.
The speed improvement is excellent.




%CONVFFT FFT-based convolution and polynomial multiplication.
%    C = CONVFFT(A, B) convolves vectors A and B.  The resulting
%    vector is length LENGTH(A)+LENGTH(B)-1.
%    If A and B are vectors of polynomial coefficients, convolving
%    them is equivalent to multiplying the two polynomials.
%
% Please contribute if you find this software useful.
% Report bugs to luigi.rosa@tiscali.it
%
%*****************************************************************

%CONV2FFT FFT-based two dimensional convolution.
%    C = CONV2FFT(A, B) performs the 2-D convolution of matrices
%    A and B.   If [ma,na] = size(A) and [mb,nb] = size(B), then
%    size(C) = [ma+mb-1,na+nb-1].
%    C = CONV2FFT(H1, H2, A) convolves A first with the vector H1 
%    along the rows and then with the vector H2 along the columns.
% 
%    C = CONV2FFT( ... ,'shape') returns a subsection of the 2-D
%    convolution with size specified by 'shape':
%      'full'  - (default) returns the full 2-D convolution,
%      'same'  - returns the central part of the convolution
%                that is the same size as A.
%      'valid' - returns only those parts of the convolution
%                that are computed without the zero-padded
%                edges. size(C) = [ma-mb+1,na-nb+1] when
%                all(size(A) >= size(B)), otherwise C is empty.
%
%
% Please contribute if you find this software useful.
% Report bugs to luigi.rosa@tiscali.it
%
%*****************************************************************
% Luigi Rosa
% Via Centrale 27
% 67042 Civita di Bagno
% L'Aquila --- ITALY 
% email  luigi.rosa@tiscali.it
% mobile +39 340 3463208 
% http://utenti.lycos.it/matlab
%*****************************************************************