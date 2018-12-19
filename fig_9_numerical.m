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
end

for n=9:length(m)
    y2=solve(1/(2*g2)*(sqrt((M(n)-1)^2+(4*g2/x)*(M(n)-1))-(M(n)+1))+1/(2*g3)*(sqrt((M(n)-1)^2+(4*g3/x)*(M(n)-1))-(M(n)+1))==Pt);
    %Solve lamda numerically
    z(n)=y2;
    
    p_1(n)=0;
    p_2(n)=1/(2*g2)*(sqrt((M(n)-1)^2+(4*g2/z(n))*(M(n)-1))-(M(n)+1));
    p_3(n)=1/(2*g3)*(sqrt((M(n)-1)^2+(4*g3/z(n))*(M(n)-1))-(M(n)+1));
end
%-------------------------------Regular WF---------------------------------
p_11=zeros(1,length(m));
p_21=zeros(1,length(m));
p_31=zeros(1,length(m));

y3=solve(2/x-1.01==Pt);   %Solve lamda numerically
p_1_wf=0;
p_2_wf=1/y3-1/g2;
p_3_wf=1/y3-1/g3;

for n=1:length(m)
    p_11(n)=p_1_wf;
    p_21(n)=p_2_wf;
    p_31(n)=p_3_wf;
end

plot(m,p_1,'linewidth',1.25)
hold on
grid on
plot(m,p_2,'r','linewidth',1.25)
plot(m,p_3,'k','linewidth',1.25)

plot(m,p_11,'--','linewidth',1.25)
plot(m,p_21,'--','linewidth',1.25)
plot(m,p_31,'--','linewidth',1.25)

xlabel('logM/log2')
ylabel('Power Allocation')
legend('g1=-20dB','g2=0dB','g3=20dB','WF')

text(6,53,'ch3')
text(6,47,'ch2')
text(6,2,'ch1')