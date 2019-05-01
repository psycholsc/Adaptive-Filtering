sigma_nr2=10;
B=[0.4];
A=[1 -1.36 0.79];
L=10000;

% var(n)
% plot(n)
% title('White Noise')
% figure
% plot(10*log10(periodogram(n)))
% title('Filtered Noise')
r_ave=zeros(1,19999);
for i = 1:10000;
    nr=sqrt(sigma_nr2)*randn(1,L);
    n=filter(B,A,nr);
    k=1:L;
    x=sin(0.2*pi*k);
    r=x+n;
    rr=xcorr(r);
    r_ave=r_ave+rr/10000;
end
r_ave=r_ave/10000;
plot(r_ave)
r_Mat=[
    r_ave(10000:10020);
    r_ave(9999:10019);
    r_ave(9998:10018);
    r_ave(9997:10017);
    r_ave(9996:10016);
    r_ave(9995:10015);
    r_ave(9994:10014);
    r_ave(9993:10013);
    r_ave(9992:10012);
    r_ave(9991:10011);
    r_ave(9990:10010);
    r_ave(9989:10009);
    r_ave(9988:10008);
    r_ave(9987:10007);
    r_ave(9986:10006);
    r_ave(9985:10005);
    r_ave(9984:10004);
    r_ave(9983:10003);
    r_ave(9982:10002);
    r_ave(9981:10001);
    r_ave(9980:10000);
    ]

%%
[A,B]=eig(r_Mat)
B(1,1)/B(21,21)
%% Begin
sigma_nr2=10;
B=[0.4];
A=[1 -1.36 0.79];
L=10000;
nr=sqrt(sigma_nr2)*randn(1,L);
n=filter(B,A,nr);
k=1:L;
x=sin(0.2*pi*k);
r=x+nr;
input=n;
ri = xcorr(input)/L;
r_ave=zeros(1,19999);
for i = 1:10000;
    nr=sqrt(sigma_nr2)*randn(1,L);
    n=filter(B,A,nr);
    r_ave=r_ave+xcorr(n)/10000;
end
r_ave=r_ave/L;
plot(r_ave);
r_Mat=[
    r_ave(10000:10020);
    r_ave(9999:10019);
    r_ave(9998:10018);
    r_ave(9997:10017);
    r_ave(9996:10016);
    r_ave(9995:10015);
    r_ave(9994:10014);
    r_ave(9993:10013);
    r_ave(9992:10012);
    r_ave(9991:10011);
    r_ave(9990:10010);
    r_ave(9989:10009);
    r_ave(9988:10008);
    r_ave(9987:10007);
    r_ave(9986:10006);
    r_ave(9985:10005);
    r_ave(9984:10004);
    r_ave(9983:10003);
    r_ave(9982:10002);
    r_ave(9981:10001);
    r_ave(9980:10000);
    ];
[A,B]=eig(r_Mat);
lambda_max=B(1,1);
%%
%   Definitions:
sigma_nr2=10;
B=[0.4];
A=[1 -1.36 0.79];
L=2000;
nr=sqrt(sigma_nr2)*randn(1,L);
n=filter(B,A,nr);
k=1:L;
x=sin(0.2*pi*k);
r=x+nr;
ensemble    = 100;                          % number of realizations within the ensemble
K           = L;                          % number of iterations
% H           = [0.32+0.21*j,-0.3+0.7*j,0.5-0.8*j,0.2+0.5*j].';
% Wo          = H;                            % unknown system
% sigma_n2    = 0.04;                         % noise power
N           = 21;                            % number of coefficients of the adaptive filter
delta       = 0.1;                          % small positive constant (used to initialize the
                                            % estimate of the inverse of the autocorrelation
                                            % matrix)
lambda      = 1.0;                         % forgetting factor


%   Initializing & Allocating memory:
W       = zeros(N,(K+1),ensemble);  % coefficient vector for each iteration and realization
MSE     = zeros(K,ensemble);        % MSE for each realization
MSEPost = zeros(K,ensemble);        % MSE a posteriori for each realization
MSEmin  = zeros(K,ensemble);        % MSE_min for each realization


%   Computing:
for l=1:ensemble,

    X       = zeros(N,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);               % desired signal
    nr      = sqrt(sigma_nr2)*randn(L,1);

    x       = filter(B,A,nr);           % Creating the input signal
                                                                % (normalized)

    % sigma_x2 = var(x);                                          % signal power = 1
    % n        = sqrt(sigma_n2/2)*(randn(K,1)+j*randn(K,1));      % complex noise
    k=1:L;
    t = [1:L]/40000;
    f0 = 1000;
    u = 10^6;
    yy = real(exp(1j*2*pi*(f0*t+1/2*u*t.*t)));
    dd=yy+transpose(nr);

    for k=1:K,

        X       =   [x(k,1)
                     X(1:(N-1),1)];              % input signal (tapped delay line)

        d(k)    =   dd(k);           % desired signal

    end

    S   =   struct('filterOrderNo',(N-1),'delta',delta,'lambda',lambda);
    [y,e,W(:,:,l),yPost,ePost]  =   RLS(d,transpose(x),S);

    MSE(:,l)    =   MSE(:,l)+(abs(e(:,1))).^2;
    MSEPost(:,l)=   MSEPost(:,l)+(abs(ePost(:,1))).^2;
    MSEmin(:,l) =   MSEmin(:,l)+(abs(n(:))).^2;

end


%   Averaging:
W_av        = sum(W,3)/ensemble;
MSE_av      = sum(MSE,2)/ensemble;
MSEPost_av  = sum(MSEPost,2)/ensemble;
MSEmin_av   = sum(MSEmin,2)/ensemble;


%   Plotting:
figure,
plot(1:K,10*log10(MSE_av),'-k');
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

figure,
plot(1:K,MSE_av,'-k');
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [Linear]');

figure,
plot(1:K,10*log10(MSEPost_av),'-k');
title('Learning Curve for MSE(a posteriori)');
xlabel('Number of iterations, k'); ylabel('MSE(a posteriori) [dB]');

figure,
plot(1:K,MSEPost_av,'-k');
title('Learning Curve for MSE(a posteriori)');
xlabel('Number of iterations, k'); ylabel('MSE(a posteriori) [Linear]');

figure,
plot(1:K,10*log10(MSEmin_av),'-k');
title('Learning Curve for MSEmin');
xlabel('Number of iterations, k'); ylabel('MSEmin [dB]');

figure,
subplot 211, plot(real(W_av(1,:))),...
title('Evolution of the 1st coefficient (real part)');
xlabel('Number of iterations, k'); ylabel('Coefficient');
subplot 212, plot(imag(W_av(1,:))),...
title('Evolution of the 1st coefficient (imaginary part)');
xlabel('Number of iterations, k'); ylabel('Coefficient');

figure
plot(e)
hold on
k=1:L;
plot(yy)
%%

figure
subplot(2,2,1);
plot(imag(y));
subplot(2,2,2);
plot(real(y));
subplot(2,2,3);
spectrogram(imag(y),256,250,256,40000);
subplot(2,2,4);
spectrogram(real(y),256,250,256,40000);
