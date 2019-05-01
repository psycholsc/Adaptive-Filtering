%% P77 e.g. 3.5
L=1500;
a=0; %0.6894, 20;0.8702, 80
R= [1,a,a^2,a^3,a^4,a^5,a^6,a^7;
    a,1,a  ,a^2,a^3,a^4,a^5,a^6;
    a^2,a,1,a  ,a^2,a^3,a^4,a^5;
    a^3,a^2,a,1,a  ,a^2,a^3,a^4;
    a^4,a^3,a^2,a,1,a  ,a^2,a^3;
    a^5,a^4,a^3,a^2,a,1,a  ,a^2;
    a^6,a^5,a^4,a^3,a^2,a,1,a  ;
    a^7,a^6,a^5,a^4,a^3,a^2,a,1];
[A,B]=eig(R);
1/B(8,8);
%%
%LMS2 Problem 1.1.1.1.2.1
%
%   'ifile.mat' - input file containing:
%      I - members of ensemble
%      K - iterations
%      sigmax - standard deviation of input
%      Wo - coefficient vector of plant
%      sigman - standard deviation of measurement noise
%      mu - convergence factor
%      b - bits in decimal part
% 
%   'ofile.mat' - output file containing:
%      ind - sample indexes 
%      MSE - mean-square error
%      MSNDW - mean-square norm of coefficient-error vector
ensemble    = 200;                          % number of realizations within the ensemble
K           = L;                          % number of iterations
H           = [0.1 0.3 0,-0.2 -0.4 -0.7 -0.4 -0.2].';
Wo          = H;                            % unknown system
% sigma_n2    = 0.04;                         % noise power
N           = 8;                            % number of coefficients of the adaptive filter
mu          = 1/B(8,8)/(N+2)/5;                           % convergence factor (step)  (0 < mu < 1)
I=ensemble;

sigma_v2=1-a^2;
sigma_n2=0.0015;
sigman=sqrt(sigma_n2);
b=16;
Lf=length(Wo);		% plant and filter length
N=Lf-1;			% plant and filter order
MSE=zeros(K,1);		% prepare to accumulate MSE*I
MSNDW=zeros(K,1);	% prepare to accumulate MSNDW*I


v=sqrt(sigma_v2)*(randn(L,1)); 
n=sqrt(sigma_n2)*(randn(L,1));       % complex noise


for i=1:I,		% ensemble
   X=zeros(Lf,1);    	% initial memory
   W=zeros(Lf,1);	% initial coefficient vector
   % x=randn(K,1)*sigmax;		% input 
   x=zeros(L,1);
    for i=1:L-1
        x(i+1)=a*x(i)+v(i);
    end
    % sigmax=sqrt(var(x))
   
   
   n=randn(K,1)*sigman;		% measurement noise 
   for k=1:K,		% iterations
      X=[x(k)		
         X(1:N)];	% new input vector
      d=Wo'*X;		% desired signal sample
      y=W'*X;		% output sample
      e=d+n(k)-y;	
      e=qround(e,b);	% error sample
      W=W+2*mu*e*X;	
      W=qround(W,b); 	% new coefficient vector
      MSE(k)=MSE(k)+e^2;	% accumulate MSE*I
      MSNDW(k)=MSNDW(k)+norm((Wo-W),2)^2; 
                        % accumulate MSNDW*I
   end
end

ind=0:(K-1);		% sample indexes
MSE=MSE/I;		% calculate MSE
MSNDW=MSNDW/I;		% calculate MSNDW
a=sum(MSE(K-99:K))/100
b=sum(MSNDW(K-99:K))/100
figure,
plot(1:K,10*log10(MSE),'-k');
title(['Learning Curve for MSE @\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('MSE [dB]');


figure,
plot(1:K,10*log10(MSNDW),'-k');
title(['Learning Curve for MSEmin @\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('MSNDW [dB]');

% save ofile ind MSE MSNDW;	% write output variables
gamma=1;
bd=16;
sigma_e2=gamma*(2^(-2*bd))/12;
sigma_w2=(2^(-2*bd))/12;