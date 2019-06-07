%% Colored Noise Gen

% White Noise

sigma=1;
ave=0;
length=8193;
WN=ave+sigma*randn(1,length);
figure;
subplot(211)
plot(WN);
xlabel('t')
ylabel('magtitude')
title('White Noise in Time Domain')
[p,w]=periodogram(WN);
subplot(212)
plot(w,10*log10(p));
xlabel('omega')
ylabel('magtitude')
title('White Noise in Freq Domain')

% Filtered (Normalized)

% ARMA(3,1)
b=[1 0.2]
a=[1 -0.2 0.1 0.3]

% ARMA(2,2)
% b=[1 0.2 -0.3]
% a=[1 -0.2 0.1]

CN=filter(b,a,WN);
CN=CN/sqrt(var(CN))*sigma;
figure;
subplot(211)
plot(CN);
xlabel('t')
ylabel('magtitude')
title('Colored Noise in Time Domain')
[p,w]=periodogram(CN);
standard=sum(p)/length;
standardp=p;
subplot(212)
plot(w,10*log10(p));
xlabel('omega')
ylabel('magtitude')
title('Colored Noise in Freq Domain')
%%
figure
zplane(b,a)
% R

ri=flip(xcorr(CN)/length);
r=xcorr(CN)/length;

%% ARMA

% AR Estimation

A_order=3;
B_order=1;
matrixHeight=A_order;
Rmatrix = zeros(matrixHeight,A_order);
for i=1:matrixHeight;
    Rmatrix(i,1:A_order)=r(-B_order+length-i+1:-B_order+length-i+A_order);
end
Rvector=r(B_order+length+1:B_order+length+matrixHeight)';
a_hat=-1*inv(Rmatrix'*Rmatrix)*Rmatrix'*Rvector

% MA Estimation
len_ahat=size(a_hat);
len_ahat=len_ahat(1);
aa=zeros(1,len_ahat);

Esti_A=[1 a_hat'];
AY=filter(Esti_A,1,CN);
AY=AY/sqrt(var(AY))*sigma;
gamma=xcorr(AY);
omega=0:1*pi/length:1*pi-1*pi/length;
% Phi=zeros(1,length);
Phi=fft(gamma)/(length*2-1);
Phi=abs(Phi(1:length));
figure
subplot(211)
plot(10*log10(Phi));
xlabel('omega')
ylabel('magtitude (dB)')
title('A(z)y(n) Power Spectrum in Freq Domain')

for w=1:length
    Aw=1;
    for i=1:len_ahat
        Aw=Aw+a_hat(i)*exp(-1j*omega(w)*i);
    end
    Phi(w)=Phi(w)/(abs(Aw)^2);
end
Phi=Phi/(sum(Phi)/length)*(standard);
subplot(212)
plot(10*log10(Phi))
axis([0 10000/pi*3.5 -60 20])
xlabel('omega')
ylabel('magtitude (dB)')
title('y(n) Power Spectrum in Freq Domain')
figure
plot(standardp)
hold on
plot(Phi)
hold off

figure
subplot(211)
plot(standardp)
title('Periodogram')
subplot(212)
plot(Phi)
title('ARMA')