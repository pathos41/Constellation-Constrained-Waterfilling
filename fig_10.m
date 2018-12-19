N=10000;
SNR=-20:0.1:25;  %in dB scale
C=zeros(length(SNR),1);
C1=zeros(length(SNR),1);
C2=zeros(length(SNR),1);
sigma=zeros(length(SNR),1);
gamma=zeros(length(SNR),1);
z=normrnd(0,1,N,1);
for a=1:length(SNR)
   sigma(a)=10^(SNR(a)/(-20));
   gamma(a)=1/sigma(a)^2;
   z1=sigma(a)*z;
   C(a)=1-mean(log(1+exp(-2*(1+z1)/sigma(a)^2))/log(2));
   C1(a)=1-log(1+exp(-gamma(a)))/log(2);
   C2(a)=1-log(1+exp(-2*gamma(a)))/log(2);
end
plot(SNR,C,'b','linewidth',1.25)
hold on
plot(SNR,C1,'--','linewidth',1.25)
plot(SNR,C2,'-.','linewidth',1.25)
grid on
xlabel('SNR[dB]')
ylabel('C[bits/s/Hz]')
legend('MC','1-log_2(1+e^{-\rho})','1-log_2(1+e^{-2\rho})')