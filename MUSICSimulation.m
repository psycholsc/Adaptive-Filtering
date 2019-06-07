clc; clear; close all;
L=1024;  % length
c = [1 -0.5];
nc = length(c) - 1;
xik=zeros(nc,1);  % init
xi=randn(L,1);  % 0 mean sigma=1

for k=1:L
    e(k)=c*[xi(k);xik];  %colored
    % update
    for i=nc:-1:2
        xik(i)=xik(i-1);
    end
    xik(1)=xi(k);
end
%===========================================
derad = pi/180;
N = 8;               % Arrays  
M = 3;               % Signals
theta = [-30 0 50];  % Directions
snr = 15;            % SNR
K = 1024;            % Quick shots

dd = 0.5;            % array distance/lambda, d=lambda/2      
d=0:dd:(N-1)*dd;     % d here is md/lambda
A=exp(-1i*2*pi*d.'*sin(theta*derad));

S=randn(M,K); 
X=A*S;
sigma2=(var(S(1,:))+var(S(2,:))+var(S(3,:)))/3/(10^(snr/10))
var(S(3,:))
WN=wgn(N,K,sigma2,'linear','complex');
e=zeros(N,K);
for i=1:N
    for j=3:K
        e(i,j)=WN(i,j)-0.7*WN(i,j-1)+0.7*WN(i,j-2);
    end
end

% X1=awgn(X,snr,'measured'); % X=AF+W
X1=X+e;
Rxx=X1*X1'/K;
%InvS=inv(Rxx); 
[EV,D]=eig(Rxx);    % eigval vec
EVA=diag(D)';
[EVA,I]=sort(EVA);  % small -> large
%EVA=fliplr(EVA);   % flip
EV=fliplr(EV(:,I)); % flip in column
for iang = 1:361    % iteration
        angle(iang)=(iang-181)/2;
        phim=derad*angle(iang);
        a=exp(-1i*2*pi*d*sin(phim)).';
        a2=exp(-1i*2*pi*d*cos(phim)).';
        L=M;    
        En=EV(:,L+1:N);
        %SP(iang)=(a'*a)/(a'*En*En'*a);
        SP(iang)=1/(a'*En*En'*a);
        b=[a';a2'];
        SP4=1./b*En*En'*b';
        SP4_1(iang)=SP4(1,1);
        SP4_2(iang)=SP4(1,2);
        SP4_3(iang)=SP4(2,1);
        SP4_4(iang)=SP4(2,2);
        
        SP2(iang)=1/(a'*inv(Rxx)*a);
        SP3(iang)=a'*Rxx*a;
end

SP=abs(SP);
%SPmax=max(SP);
%SP=10*log10(SP/SPmax);  % normalized
SP=10*log10(SP);
h=plot(angle,SP);
set(h,'Linewidth',0.5);
xlabel('DOA/(degree)');
ylabel('Spectrum/(dB)');
%axis([-100 100 -40 60]);
set(gca, 'XTick',[-100:20:100]);
grid on
hold on;

%%
SP2=abs(SP2);
%SPmax=max(SP);
%SP=10*log10(SP/SPmax);  % normalized
SP2=10*log10(SP2);
h2=plot(angle,SP2);
set(h2,'Linewidth',0.5);
xlabel('DOA/(degree)');
ylabel('Spectrum/(dB)');
%axis([-100 100 -40 60]);
set(gca, 'XTick',[-100:20:100]);
grid on
hold on;

SP3=abs(SP3);
%SPmax=max(SP);
%SP=10*log10(SP/SPmax);  % normalized
SP3=10*log10(SP3);
h3=plot(angle,SP3);
set(h3,'Linewidth',0.5);
xlabel('DOA/(degree)');
ylabel('Spectrum/(dB)');
%axis([-100 100 -40 60]);
set(gca, 'XTick',[-100:20:100]);
grid on
legend('MUSIC','ML','BF')

%=======================
figure
[y1,f1] = periodogram(e(1,:));
p1 = 10*log10(y1);
plot(f1,p1)
xlabel('freq')
ylabel('magtitude (dB)')
title('periodogram of colored noise')