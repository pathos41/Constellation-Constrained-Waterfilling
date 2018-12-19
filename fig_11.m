%----------------------------Initial parameters----------------------------
d_2=[0:1;
    -1:0]'*2;  % Distance matrix
d_4=[0:3;
    -1:2;
    -2:1;
    -3:0]'*2/sqrt(5);

N=10000;
SNR=-20:0.1:30;
sigma=10.^(SNR/(-20));  % Standard deviation of noise
gamma=1./sigma.^2;

C_2=zeros(length(SNR),1);
C_4=zeros(length(SNR),1);

f1=zeros(length(d_2),1);
f2=zeros(length(d_4),1);

b1=zeros(1,length(d_2));
b2=zeros(1,length(d_4));

x1=zeros(N,1);
x2=zeros(N,1);

C1_2=zeros(length(SNR),1);
C1_4=zeros(length(SNR),1);

C2_2=zeros(length(SNR),1);
C2_4=zeros(length(SNR),1);
%----------------------------------MC--------------------------------------
z=normrnd(0,1,N,1);  % Gaussian random variable
%----------------------------------M=2-------------------------------------
for a=1:length(SNR)
    z1=sigma(a)*z;  % Gaussian random variable
    for m=1:N  % Number of pages
        for j=1:length(d_2)  % Number of columns
            for i=1:length(d_2)  % Number of rows
                f1(i)=exp(-d_2(i,j)^2/(2*sigma(a)^2)-z1(m)*d_2(i,j)/sigma(a)^2);
            end
            b1(j)=log(sum(f1))/log(2);  % First sum over i
        end
        x1(m)=sum(b1);  % Second sum over j
    end
    C_2(a)=1-(1/2)*mean(x1);  % Capacity of 2-PAM
end
%----------------------------------M=4-------------------------------------
for a=1:length(SNR)
    z2=sigma(a)*z;  % Gaussian random variable
    for m=1:N  % Number of pages
        for j=1:length(d_4)  % Number of columns
            for i=1:length(d_4)  % Number of rows
                f2(i)=exp(-d_4(i,j)^2/(2*sigma(a)^2)-z2(m)*d_4(i,j)/sigma(a)^2);
            end
            b2(j)=log(sum(f2))/log(2);  % First sum over i
        end
        x2(m)=sum(b2);  % Second sum over j
    end
    C_4(a)=2-(1/4)*mean(x2);  % Capacity of 4-PAM
end
%----------------------------Sphere Packing--------------------------------
for a=1:length(SNR)
    z=normrnd(0,sigma(a),N,1);
    C1_2(a)=0.5*log((1+gamma(a))/(1+gamma(a)/2^2))/log(2);
    C1_4(a)=0.5*log((1+gamma(a))/(1+gamma(a)/4^2))/log(2);
end
%------------------------------Analytical----------------------------------
for a=1:length(SNR)
    C2_2(a)=1-log(1+exp(-2*gamma(a)))/log(2);
    C2_4(a)=2-log(1+exp((-2/5)*gamma(a)))/log(2);
end

plot(SNR,C_2,'linewidth',1.25)
hold on
plot(SNR,C1_2,'r','linewidth',1.25)
plot(SNR,C2_2,'k','linewidth',1.25)

plot(SNR,C_4,'linewidth',1.25)
plot(SNR,C1_4,'r','linewidth',1.25)
plot(SNR,C2_4,'k','linewidth',1.25)

grid on

xlabel('SNR[dB]')
ylabel('C[bits/s/Hz]')

text(15,0.7,'2-PAM')
text(20,1.7,'4-PAM')

legend('MC','Sphere-Packing','High-SNR Approx')