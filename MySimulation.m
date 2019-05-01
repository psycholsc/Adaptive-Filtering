%% P77 e.g. 3.5

%% Question (a)+(b)

L           =   1500;               % Iterations and Siganal length
a           =   0;                  %a = 0, 1; a = 0.6894, 20; a = 0.8702, 80
sigma_v2    =   1-a^2;              % Make sure colored noise power = 1
sigma_n2    =   0.0001;             % measure noise
R           =   [1,a,a^2,a^3,a^4,a^5,a^6,a^7;
                a,1,a  ,a^2,a^3,a^4,a^5,a^6;
                a^2,a,1,a  ,a^2,a^3,a^4,a^5;
                a^3,a^2,a,1,a  ,a^2,a^3,a^4;
                a^4,a^3,a^2,a,1,a  ,a^2,a^3;
                a^5,a^4,a^3,a^2,a,1,a  ,a^2;
                a^6,a^5,a^4,a^3,a^2,a,1,a  ;
                a^7,a^6,a^5,a^4,a^3,a^2,a,1];
[A,B]       =   eig(R);


v=sqrt(sigma_v2)*(randn(L,1));      % White Noise
n=sqrt(sigma_n2)*(randn(L,1));

plot(10*log10(periodogram(v)),'-k');
title('Periodogram of White Noise')
x=zeros(L,1);
for i=1:L-1
    x(i+1)=a*x(i)+v(i);             % Generate Colored Noise
end
figure
plot(10*log10(periodogram(x)),'-k');
title(['Periodogram of Colored Noise @a=' num2str(a)])

%   Definitions:
ensemble    = 100;                          % number of realizations within the ensemble
K           = L;                            % number of iterations
H           = [0.1 0.3 0,-0.2 -0.4 -0.7 -0.4 -0.2].';
Wo          = H;                            % unknown system
N           = 8;                            % number of coefficients of the adaptive filter
mu          = 1/B(8,8)/(N+2)                           % convergence factor (step)  (0 < mu < 1)

%   Initializing & Allocating memory:
W       = ones(N,(K+1),ensemble);   % coefficient vector for each iteration and realization; w(0) = [1 1 1 1].'
MSE     = zeros(K,ensemble);        % MSE for each realization
MSEmin  = zeros(K,ensemble);        % MSE_min for each realization
M=0;
%   Computing:
for l=1:ensemble,

    X       = zeros(N,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);               % desired signal
    n=sqrt(sigma_n2)*(randn(L,1));
    v=sqrt(sigma_v2)*(randn(L,1)); 
    x=zeros(L,1);
    for i=1:L-1
        x(i+1)=a*x(i)+v(i);
    end
    % nw=sqrt(sigma_w2)*randn(1,K);
    
    for k=1:K,
        X       =   [x(k,1)
                     X(1:(N-1),1)];              % input signal (tapped delay line)
        % Wo=lambda*Wo+nw(k);                      % For Non-stationary
        d(k)    =   (Wo'*X(:,1))+n(k);           % desired signal
    end
    S   =   struct('step',mu,'filterOrderNo',(N-1),'initialCoefficients',W(:,1,l));
    [y,e,W(:,:,l)]  =   LMS(d,transpose(x),S);
    MSE(:,l)    =   MSE(:,l)+(abs(e(:,1))).^2;
    MSEmin(:,l) =   MSEmin(:,l)+(abs(n(:))).^2;
    M=M+mu*8*var(x)/(1-mu*8*var(x));

end
%   Averaging:
W_av = sum(W,3)/ensemble;
MSE_av = sum(MSE,2)/ensemble;
MSEmin_av = sum(MSEmin,2)/ensemble;

Mreal=(sum(MSE_av(K-99:K))/100-mean(MSEmin_av))/mean(MSEmin_av)

M=M/ensemble

%   Plotting:
figure,
plot(1:K,10*log10(MSE_av),'-k');
title(['Learning Curve for MSE @\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

figure,
plot(1:K,10*log10(MSEmin_av),'-k');
title(['Learning Curve for MSEmin @\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('MSEmin [dB]');

figure,
subplot 211, plot(real(W_av(1,:))),...
title(['Evolution of the 1st coefficient (real part)@\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('Coefficient');
subplot 212, plot(imag(W_av(1,:))),...
title(['Evolution of the 1st coefficient (imaginary part)@\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('Coefficient');



lambda=0.99; % Nonstationary
sigma_w2=0.0015;


%%
%








%%
ensemble    = 150;                          % number of realizations within the ensemble
K           = L;                          % number of iterations
H           = [0.1 0.3 0,-0.2 -0.4 -0.7 -0.4 -0.2].';
Wo          = H;                            % unknown system
% sigma_n2    = 0.04;                         % noise power
N           = 8;                            % number of coefficients of the adaptive filter
mu          = 1/B(8,8)/(N+2)/5;                           % convergence factor (step)  (0 < mu < 1)


%   Initializing & Allocating memory:
W       = ones(N,(K+1),ensemble);   % coefficient vector for each iteration and realization; w(0) = [1 1 1 1].'
MSE     = zeros(K,ensemble);        % MSE for each realization
MSEmin  = zeros(K,ensemble);        % MSE_min for each realization
M=0;

%   Computing:
for l=1:ensemble,

    X       = zeros(N,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);               % desired signal

%     x        = (sign(randn(K,1)) + j*sign(randn(K,1)))./sqrt(2); % Creating the input signal (normalized)
%     var(x) = var(x);                                           % signal power = 1
%     n        = sqrt(sigma_n2/2)*(randn(K,1)+j*randn(K,1));       % complex noise
    n=sqrt(sigma_n2)*(randn(L,1));
    v=sqrt(sigma_v2)*(randn(L,1)); 
    x=zeros(L,1);
    for i=1:L-1
        x(i+1)=a*x(i)+v(i);
    end
    
    for k=1:K,

        X       =   [x(k,1)
                     X(1:(N-1),1)];              % input signal (tapped delay line)

        d(k)    =   (Wo'*X(:,1))+n(k);           % desired signal

    end

    S   =   struct('step',mu,'filterOrderNo',(N-1),'initialCoefficients',W(:,1,l));
    [y,e,W(:,:,l)]  =   LMS(d,transpose(x),S);

    MSE(:,l)    =   MSE(:,l)+(abs(e(:,1))).^2;
    MSEmin(:,l) =   MSEmin(:,l)+(abs(n(:))).^2;
    M=M+mu*8*var(x)/(1-mu*8*var(x));

end


%   Averaging:
W_av = sum(W,3)/ensemble;
MSE_av = sum(MSE,2)/ensemble;
MSEmin_av = sum(MSEmin,2)/ensemble;
Mreal=(mean(MSEmin_av)-sigma_n2)/sigma_n2

%   Plotting:
figure,
plot(1:K,10*log10(MSE_av),'-k');
title(['Learning Curve for MSE @\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

figure,
plot(1:K,10*log10(MSEmin_av),'-k');
title(['Learning Curve for MSEmin @\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('MSEmin [dB]');

figure,
subplot 211, plot(real(W_av(1,:))),...
title(['Evolution of the 1st coefficient (real part)@\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('Coefficient');
subplot 212, plot(imag(W_av(1,:))),...
title(['Evolution of the 1st coefficient (imaginary part)@\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('Coefficient');
M=M/ensemble
%%
ensemble    = 100;                          % number of realizations within the ensemble
K           = L;                          % number of iterations
H           = [0.1 0.3 0,-0.2 -0.4 -0.7 -0.4 -0.2].';
Wo          = H;                            % unknown system
% sigma_n2    = 0.04;                         % noise power
N           = 8;                            % number of coefficients of the adaptive filter
mu          = 1/B(8,8)/(N+2)/10;                           % convergence factor (step)  (0 < mu < 1)


%   Initializing & Allocating memory:
W       = ones(N,(K+1),ensemble);   % coefficient vector for each iteration and realization; w(0) = [1 1 1 1].'
MSE     = zeros(K,ensemble);        % MSE for each realization
MSEmin  = zeros(K,ensemble);        % MSE_min for each realization
M=0;

%   Computing:
for l=1:ensemble,

    X       = zeros(N,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);               % desired signal

%     x        = (sign(randn(K,1)) + j*sign(randn(K,1)))./sqrt(2); % Creating the input signal (normalized)
%     var(x) = var(x);                                           % signal power = 1
%     n        = sqrt(sigma_n2/2)*(randn(K,1)+j*randn(K,1));       % complex noise
    n=sqrt(sigma_n2)*(randn(L,1));
    v=sqrt(sigma_v2)*(randn(L,1)); 
    x=zeros(L,1);
    for i=1:L-1
        x(i+1)=a*x(i)+v(i);
    end
    
    for k=1:K,

        X       =   [x(k,1)
                     X(1:(N-1),1)];              % input signal (tapped delay line)

        d(k)    =   (Wo'*X(:,1))+n(k);           % desired signal

    end

    S   =   struct('step',mu,'filterOrderNo',(N-1),'initialCoefficients',W(:,1,l));
    [y,e,W(:,:,l)]  =   LMS(d,transpose(x),S);

    MSE(:,l)    =   MSE(:,l)+(abs(e(:,1))).^2;
    MSEmin(:,l) =   MSEmin(:,l)+(abs(n(:))).^2;
    M=M+mu*8*var(x)/(1-mu*8*var(x));

end


%   Averaging:
W_av = sum(W,3)/ensemble;
MSE_av = sum(MSE,2)/ensemble;
MSEmin_av = sum(MSEmin,2)/ensemble;
Mreal=(mean(MSEmin_av)-sigma_n2)/sigma_n2

%   Plotting:
figure,
plot(1:K,10*log10(MSE_av),'-k');
title(['Learning Curve for MSE @\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

figure,
plot(1:K,10*log10(MSEmin_av),'-k');
title(['Learning Curve for MSEmin @\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('MSEmin [dB]');

figure,
subplot 211, plot(real(W_av(1,:))),...
title(['Evolution of the 1st coefficient (real part)@\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('Coefficient');
subplot 212, plot(imag(W_av(1,:))),...
title(['Evolution of the 1st coefficient (imaginary part)@\mu=' num2str(mu)]);
xlabel('Number of iterations, k'); ylabel('Coefficient');
M=M/ensemble