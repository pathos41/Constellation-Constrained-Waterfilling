g1=1/100;
g2=1;
g3=100;  %Channel power gains

m=1:25;  %logM/log2
M=2.^m;  %M-symbol constellation
Pt=100;  %Total power

c_1=zeros(1,length(m));
c_2=zeros(1,length(m));
c_3=zeros(1,length(m));

for n=1:length(m)
    cvx_begin  %Get optimal power allocation
    cvx_precision low
    variable p1
    variable p2
    variable p3
    
    a=log(1+(M(n)-1)*(1-M(n)*inv_pos(g1*p1+M(n))));  %Constellation constrained capacity
    b=log(1+(M(n)-1)*(1-M(n)*inv_pos(g2*p2+M(n))));
    c=log(1+(M(n)-1)*(1-M(n)*inv_pos(g3*p3+M(n))));
    
    minimize(-a-b-c);
    
    subject to
    p1+p2+p3<=Pt  %Power constraints
    -p1<=0
    -p2<=0
    -p3<=0
    cvx_end
    c_1(n)=a/log(2);
    c_2(n)=b/log(2);
    c_3(n)=c/log(2);
end

plot(m,c_1,'linewidth',1.25)
hold on
grid on
plot(m,c_2,'r','linewidth',1.25)
plot(m,c_3,'k','linewidth',1.25)

xlim([0 25])
ylim([0 14])

xlabel('log_2M')
ylabel('Constellation capacity [bit/s/Hz]')
legend('g_1=-20dB','g_2=0dB','g_3=20dB')