%----------------------------Distance Matrices-----------------------------
d_2=[0:1;
    -1:0]'*2;  % Distance matrix
d_4=[0:3;
    -1:2;
    -2:1;
    -3:0]'*2/sqrt(5);
d_8=[0:7;
    -1:6;
    -2:5;
    -3:4;
    -4:3;
    -5:2;
    -6:1;
    -7:0]'*2/sqrt(21);
d_16=[0:15;
    -1:14;
    -2:13;
    -3:12;
    -4:11;
    -5:10;
    -6:9;
    -7:8;
    -8:7;
    -9:6;
    -10:5;
    -11:4;
    -12:3;
    -13:2;
    -14:1;
    -15:0]'*2/sqrt(85);
d_32=[0:31;
    -1:30;
    -2:29;
    -3:28;
    -4:27;
    -5:26;
    -6:25;
    -7:24;
    -8:23;
    -9:22;
    -10:21;
    -11:20;
    -12:19;
    -13:18;
    -14:17;
    -15:16;
    -16:15;
    -17:14;
    -18:13;
    -19:12;
    -20:11;
    -21:10;
    -22:9;
    -23:8;
    -24:7;
    -25:6;
    -26:5;
    -27:4;
    -28:3;
    -29:2;
    -30:1;
    -31:0]'*2/sqrt(341);
d_64=[0:63;
    -1:62;
    -2:61;
    -3:60;
    -4:59;
    -5:58;
    -6:57;
    -7:56;
    -8:55;
    -9:54;
    -10:53;
    -11:52;
    -12:51;
    -13:50;
    -14:49;
    -15:48;
    -16:47;
    -17:46;
    -18:45;
    -19:44;
    -20:43;
    -21:42;
    -22:41;
    -23:40;
    -24:39;
    -25:38;
    -26:37;
    -27:36;
    -28:35;
    -29:34;
    -30:33;
    -31:32;
    -32:31;
    -33:30;
    -34:29;
    -35:28;
    -36:27;
    -37:26;
    -38:25;
    -39:24;
    -40:23;
    -41:22;
    -42:21;
    -43:20;
    -44:19;
    -45:18;
    -46:17;
    -47:16;
    -48:15;
    -49:14;
    -50:13;
    -51:12;
    -52:11;
    -53:10;
    -54:9;
    -55:8;
    -56:7;
    -57:6;
    -58:5;
    -59:4;
    -60:3;
    -61:2;
    -62:1;
    -63:0]'*2/sqrt(1365);
%----------------------------Initial parameters----------------------------
N=10000;
SNR=-20:60;
sigma=10.^(SNR/(-20));  % Standard deviation of noise
gamma=1./sigma.^2;

C_2=zeros(length(SNR),1);
C_4=zeros(length(SNR),1);
C_8=zeros(length(SNR),1);
C_16=zeros(length(SNR),1);
C_32=zeros(length(SNR),1);
C_64=zeros(length(SNR),1);

f1=zeros(length(d_2),1);
f2=zeros(length(d_4),1);
f3=zeros(length(d_8),1);
f4=zeros(length(d_16),1);
f5=zeros(length(d_32),1);
f6=zeros(length(d_64),1);

b1=zeros(1,length(d_2));
b2=zeros(1,length(d_4));
b3=zeros(1,length(d_8));
b4=zeros(1,length(d_16));
b5=zeros(1,length(d_32));
b6=zeros(1,length(d_64));

x1=zeros(N,1);
x2=zeros(N,1);
x3=zeros(N,1);
x4=zeros(N,1);
x5=zeros(N,1);
x6=zeros(N,1);

C=zeros(length(SNR),1);
C1_2=zeros(length(SNR),1);
C1_4=zeros(length(SNR),1);
C1_8=zeros(length(SNR),1);
C1_16=zeros(length(SNR),1);
C1_32=zeros(length(SNR),1);
C1_64=zeros(length(SNR),1);

for a=1:length(SNR)
    z=normrnd(0,sigma(a),N,1);
    C(a)=0.5*log(1+gamma(a))/log(2);
    C1_2(a)=0.5*log((1+gamma(a))/(1+gamma(a)/2^2))/log(2);
    C1_4(a)=0.5*log((1+gamma(a))/(1+gamma(a)/4^2))/log(2);
    C1_8(a)=0.5*log((1+gamma(a))/(1+gamma(a)/8^2))/log(2);
    C1_16(a)=0.5*log((1+gamma(a))/(1+gamma(a)/16^2))/log(2);
    C1_32(a)=0.5*log((1+gamma(a))/(1+gamma(a)/32^2))/log(2);
    C1_64(a)=0.5*log((1+gamma(a))/(1+gamma(a)/64^2))/log(2);
end
%----------------------------------M=2-------------------------------------
for a=1:length(SNR)
    z=normrnd(0,sigma(a),N,1);  % Gaussian random variable
    for m=1:N  % Number of pages
        for j=1:length(d_2)  % Number of columns
            for i=1:length(d_2)  % Number of rows
                f1(i)=exp(-d_2(i,j)^2/(2*sigma(a)^2)-z(m)*d_2(i,j)/sigma(a)^2);
            end
            b1(j)=log(sum(f1))/log(2);  % First sum over i
        end
        x1(m)=sum(b1);  % Second sum over j
    end
    C_2(a)=1-(1/2)*mean(x1);  % Capacity of 2-PAM
end
%----------------------------------M=4-------------------------------------
for a=1:length(SNR)
    z=normrnd(0,sigma(a),N,1);  % Gaussian random variable
    for m=1:N  % Number of pages
        for j=1:length(d_4)  % Number of columns
            for i=1:length(d_4)  % Number of rows
                f2(i)=exp(-d_4(i,j)^2/(2*sigma(a)^2)-z(m)*d_4(i,j)/sigma(a)^2);
            end
            b2(j)=log(sum(f2))/log(2);  % First sum over i
        end
        x2(m)=sum(b2);  % Second sum over j
    end
    C_4(a)=2-(1/4)*mean(x2);  % Capacity of 4-PAM
end
%----------------------------------M=8-------------------------------------
for a=1:length(SNR)
    z=normrnd(0,sigma(a),N,1);  % Gaussian random variable
    for m=1:N  % Number of pages
        for j=1:length(d_8)  % Number of columns
            for i=1:length(d_8)  % Number of rows
                f3(i)=exp(-d_8(i,j)^2/(2*sigma(a)^2)-z(m)*d_8(i,j)/sigma(a)^2);
            end
            b3(j)=log(sum(f3))/log(2);  % First sum over i
        end
        x3(m)=sum(b3);  % Second sum over j
    end
    C_8(a)=3-(1/8)*mean(x3);  % Capacity of 8-PAM
end
%----------------------------------M=16------------------------------------
for a=1:length(SNR)
    z=normrnd(0,sigma(a),N,1);  % Gaussian random variable
    for m=1:N  % Number of pages
        for j=1:length(d_16)  % Number of columns
            for i=1:length(d_16)  % Number of rows
                f4(i)=exp(-d_16(i,j)^2/(2*sigma(a)^2)-z(m)*d_16(i,j)/sigma(a)^2);
            end
            b4(j)=log(sum(f4))/log(2);  % First sum over i
        end
        x4(m)=sum(b4);  % Second sum over j
    end
    C_16(a)=4-(1/16)*mean(x4);  % Capacity of 16-PAM
end
%----------------------------------M=32------------------------------------
for a=1:length(SNR)
    z=normrnd(0,sigma(a),N,1);  % Gaussian random variable
    for m=1:N  % Number of pages
        for j=1:length(d_32)  % Number of columns
            for i=1:length(d_32)  % Number of rows
                f5(i)=exp(-d_32(i,j)^2/(2*sigma(a)^2)-z(m)*d_32(i,j)/sigma(a)^2);
            end
            b5(j)=log(sum(f5))/log(2);  % First sum over i
        end
        x5(m)=sum(b5);  % Second sum over j
    end
    C_32(a)=5-(1/32)*mean(x5);  % Capacity of 32-PAM
end
%----------------------------------M=64------------------------------------
for a=1:length(SNR)
    z=normrnd(0,sigma(a),N,1);  % Gaussian random variable
    for m=1:N  % Number of pages
        for j=1:length(d_64)  % Number of columns
            for i=1:length(d_64)  % Number of rows
                f6(i)=exp(-d_64(i,j)^2/(2*sigma(a)^2)-z(m)*d_64(i,j)/sigma(a)^2);
            end
            b6(j)=log(sum(f6))/log(2);  % First sum over i
        end
        x6(m)=sum(b6);  % Second sum over j
    end
    C_64(a)=6-(1/64)*mean(x6);  % Capacity of 64-PAM
end

plot(SNR,C,'--','linewidth',1.25)
hold on
grid on
plot(SNR,C_2,'r','linewidth',1.25)
plot(SNR,C1_2,'k','linewidth',1.25)

plot(SNR,C_4,'r','linewidth',1.25)
plot(SNR,C_8,'r','linewidth',1.25)
plot(SNR,C_16,'r','linewidth',1.25)
plot(SNR,C_32,'r','linewidth',1.25)
plot(SNR,C_64,'r','linewidth',1.25)


plot(SNR,C1_4,'k','linewidth',1.25)
plot(SNR,C1_8,'k','linewidth',1.25)
plot(SNR,C1_16,'k','linewidth',1.25)
plot(SNR,C1_32,'k','linewidth',1.25)
plot(SNR,C1_64,'k','linewidth',1.25)

xlabel('SNR[dB]')
ylabel('C[bits/s/Hz]')

text(40,1.2,'2-PAM')
text(40,2.2,'4-PAM')
text(40,3.2,'8-PAM')
text(40,4.2,'16-PAM')
text(40,5.2,'32-PAM')
text(40,6.2,'64-PAM')

legend('C','C_M','C_C')