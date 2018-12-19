g1=1/100;
g2=1;
g3=100;  %Channel power gains

m=1:25;  %logM/log2
M=2.^m;  %M-symbol constellation
Pt=100;  %Total power

p_1=zeros(1,length(m));
p_2=zeros(1,length(m));
p_3=zeros(1,length(m));

for n=1:length(m)
    cvx_begin  %Get optimal power allocation
    cvx_solver SeDuMi
    cvx_precision low
    variable p1
    variable p2
    variable p3
    
    a=log(1+(M(n)-1)*(1-M(n)*inv_pos(g1*p1+M(n))));  %Constellation constrained capacity
    b=log(1+(M(n)-1)*(1-M(n)*inv_pos(g2*p2+M(n))));
    c=log(1+(M(n)-1)*(1-M(n)*inv_pos(g3*p3+M(n))));
    
    minimize(-a-b-c)
    
    subject to
    p1+p2+p3==Pt;  %Power constraints
    -p1<=0;
    -p2<=0;
    -p3<=0;
    cvx_end
    
    p_1(n)=p1;
    p_2(n)=p2;
    p_3(n)=p3;
end
cs1=spline(m,p_1);
xx1=1:0.1:25;
plot(xx1,ppval(cs1,xx1),'linewidth',1.25)
hold on
grid on
cs2=spline(m,p_2);
xx2=1:0.1:25;
plot(xx2,ppval(cs2,xx2),'--','linewidth',1.25)
cs3=spline(m,p_3);
xx3=1:0.1:25;
plot(xx3,ppval(cs3,xx3),'-.','linewidth',1.25)
%plot(m,p_1,'linewidth',1.25)
%plot(m,p_2,'r','linewidth',1.25)
%plot(m,p_3,'k','linewidth',1.25)

xlabel('log_2M')
ylabel('Power Allocation')
legend('g_1=-20dB','g_2=0dB','g_3=20dB')