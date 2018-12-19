syms x

g1=1/100;
g2=1;
g3=100;  %Channel power gains

m=1:25;  %logM/log2
M=2.^m;  %M-symbol constellation
Pt=100;  %Total power

p_1=zeros(1,length(m));
p_2=zeros(1,length(m));
p_3=zeros(1,length(m));

c1=zeros(1,length(m));
c2=zeros(1,length(m));
c3=zeros(1,length(m));

y1=zeros(4,1);
z=zeros(1,length(m));
%----------------------------Constellation WF------------------------------
for n=1:8
    y1=solve(1/(2*g1)*(sqrt((M(n)-1)^2+(4*g1/x)*(M(n)-1))-(M(n)+1))+1/(2*g2)*(sqrt((M(n)-1)^2+(4*g2/x)*(M(n)-1))-(M(n)+1))+1/(2*g3)*(sqrt((M(n)-1)^2+(4*g3/x)*(M(n)-1))-(M(n)+1))==Pt);
    %Solve lamda numerically
    z(n)=y1(4);
    
    p_1(n)=1/(2*g1)*(sqrt((M(n)-1)^2+(4*g1/z(n))*(M(n)-1))-(M(n)+1));
    p_2(n)=1/(2*g2)*(sqrt((M(n)-1)^2+(4*g2/z(n))*(M(n)-1))-(M(n)+1));
    p_3(n)=1/(2*g3)*(sqrt((M(n)-1)^2+(4*g3/z(n))*(M(n)-1))-(M(n)+1));
    
    c1(n)=log((1+p_1(n)*g1)/(1+p_1(n)*g1/M(n)))/log(2);
    c2(n)=log((1+p_2(n)*g2)/(1+p_2(n)*g2/M(n)))/log(2);
    c3(n)=log((1+p_3(n)*g3)/(1+p_3(n)*g3/M(n)))/log(2);
end

for n=9:length(m)
    y2=solve(1/(2*g2)*(sqrt((M(n)-1)^2+(4*g2/x)*(M(n)-1))-(M(n)+1))+1/(2*g3)*(sqrt((M(n)-1)^2+(4*g3/x)*(M(n)-1))-(M(n)+1))==Pt);
    %Solve lamda numerically
    z(n)=y2;
    
    p_1(n)=0;
    p_2(n)=1/(2*g2)*(sqrt((M(n)-1)^2+(4*g2/z(n))*(M(n)-1))-(M(n)+1));
    p_3(n)=1/(2*g3)*(sqrt((M(n)-1)^2+(4*g3/z(n))*(M(n)-1))-(M(n)+1));
    
    c1(n)=log((1+p_1(n)*g1)/(1+p_1(n)*g1/M(n)))/log(2);
    c2(n)=log((1+p_2(n)*g2)/(1+p_2(n)*g2/M(n)))/log(2);
    c3(n)=log((1+p_3(n)*g3)/(1+p_3(n)*g3/M(n)))/log(2);
end
%-------------------------------Regular WF---------------------------------
c_11=zeros(1,length(m));
c_21=zeros(1,length(m));
c_31=zeros(1,length(m));

y3=solve(2/x-1.01==Pt);   %Solve lamda numerically
p_1_wf=0;
p_2_wf=1/y3-1/g2;
p_3_wf=1/y3-1/g3;

c11=log(1+p_1_wf*g1)/log(2);
c21=log(1+p_2_wf*g2)/log(2);
c31=log(1+p_3_wf*g3)/log(2);

for n=1:length(m)
    c_11(n)=c11;
    c_21(n)=c21;
    c_31(n)=c31;
end

plot(m,c1,'linewidth',1.25)
hold on
grid on
plot(m,c2,'r','linewidth',1.25)
plot(m,c3,'k','linewidth',1.25)

plot(m,c_11,'--','linewidth',1.25)
plot(m,c_21,'--','linewidth',1.25)
plot(m,c_31,'--','linewidth',1.25)

xlabel('log_2M')
ylabel('Constellation capacity [bit/s/Hz]')
legend('g_1=-20dB','g_2=0dB','g_3=20dB','WF')

text(20,12.7,'ch3')
text(20,6,'ch2')
text(20,0.5,'ch1')