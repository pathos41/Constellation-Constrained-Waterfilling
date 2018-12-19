g1=1/100;
g2=1;
g3=100;  %Channel power gains

m=1:25;  %logM/log2
M=2.^m;  %M-symbol constellation
Pt=100;  %Total power

p_1=zeros(1,length(m));
p_2=zeros(1,length(m));
p_3=zeros(1,length(m));

p_11=zeros(1,length(m));
p_21=zeros(1,length(m));
p_31=zeros(1,length(m));

a=zeros(1,length(m));
b=zeros(1,length(m));
max1=zeros(1,length(m));
tol=eps;   %Tolerance
%-----------------------------Constellation WF-----------------------------
for n=1:length(m)
    a(n)=0;   %Original left side
    b(n)=g3*(1-1/M(n));   %Original right side
    max1(n)=-1+ceil((log(b(n)-a(n))-log(tol))/log(2));   %Number of iterations
    for k=1:max1(n)+1
        lamda=zeros(1,length(k));
        
        p_1a=zeros(1,length(k));
        p_2a=zeros(1,length(k));
        p_3a=zeros(1,length(k));
        
        p_1b=zeros(1,length(k));
        p_2b=zeros(1,length(k));
        p_3b=zeros(1,length(k));
        
        p_1l=zeros(1,length(k));
        p_2l=zeros(1,length(k));
        p_3l=zeros(1,length(k));
        
        fa=zeros(1,length(k));
        fb=zeros(1,length(k));
        fl=zeros(1,length(k));
        
        lamda(k)=(a(n)+b(n))/2;   %Bisection
        p_1a(k)=constellation(g1,M(n),a(n));   %Power allocation at a
        p_2a(k)=constellation(g2,M(n),a(n));
        p_3a(k)=constellation(g3,M(n),a(n));
        
        p_1b(k)=constellation(g1,M(n),b(n));   %Power allocation at b
        p_2b(k)=constellation(g2,M(n),b(n));
        p_3b(k)=constellation(g3,M(n),b(n));
        
        p_1l(k)=constellation(g1,M(n),lamda(k));   %Power allocation at lamda
        p_2l(k)=constellation(g2,M(n),lamda(k));
        p_3l(k)=constellation(g3,M(n),lamda(k));
        
        fa(k)=p_1a(k)+p_2a(k)+p_3a(k)-Pt;   %Total power is 100
        fb(k)=p_1b(k)+p_2b(k)+p_3b(k)-Pt;
        fl(k)=p_1l(k)+p_2l(k)+p_3l(k)-Pt;
        
        if fl(k)==0   %Find the optimal point
            p_1(n)=p_1l(k);
            p_2(n)=p_2l(k);
            p_3(n)=p_3l(k);
            break
            
        else if fb(k)*fl(k)>0
                b(n)=lamda(k);
                
            else
                a(n)=lamda(k);
            end
        end
        if b(n)-a(n)<tol
            p_1(n)=p_1l(k);
            p_2(n)=p_2l(k);
            p_3(n)=p_3l(k);
            break
        end
    end
end
%--------------------------------Regular WF--------------------------------
c=0;   %Original left side
d=g3;   %Original right side
max2=-1+ceil((log(d-c)-log(tol))/log(2));  %Number of iterations

for kk=1:max2+1
        lamda_new=zeros(1,length(kk));
        
        p_1c=zeros(1,length(kk));
        p_2c=zeros(1,length(kk));
        p_3c=zeros(1,length(kk));
        
        p_1d=zeros(1,length(kk));
        p_2d=zeros(1,length(kk));
        p_3d=zeros(1,length(kk));
        
        p_1l_new=zeros(1,length(kk));
        p_2l_new=zeros(1,length(kk));
        p_3l_new=zeros(1,length(kk));
        
        fc=zeros(1,length(kk));
        fd=zeros(1,length(kk));
        fl_new=zeros(1,length(kk));
        
        lamda_new(kk)=(c+d)/2;   %Bisection
        p_1c(kk)=regular_wf(g1,c);   %Power allocation at c
        p_2c(kk)=regular_wf(g2,c);
        p_3c(kk)=regular_wf(g3,c);
        
        p_1d(kk)=regular_wf(g1,d);   %Power allocation at d
        p_2d(kk)=regular_wf(g2,d);
        p_3d(kk)=regular_wf(g3,d);
        
        p_1l_new(kk)=regular_wf(g1,lamda_new(kk));   %Power allocation at lamda
        p_2l_new(kk)=regular_wf(g2,lamda_new(kk));
        p_3l_new(kk)=regular_wf(g3,lamda_new(kk));
        
        fc(kk)=p_1c(kk)+p_2c(kk)+p_3c(kk)-Pt;   %Total power is 100
        fd(kk)=p_1d(kk)+p_2d(kk)+p_3d(kk)-Pt;
        fl_new(kk)=p_1l_new(kk)+p_2l_new(kk)+p_3l_new(kk)-Pt;
        
        if fl_new(kk)==0   %Find the optimal point
            p_1_wf=p_1l_new(kk);
            p_2_wf=p_2l_new(kk);
            p_3_wf=p_3l_new(kk);
            break
            
        else if fd(kk)*fl_new(kk)>0
                d=lamda_new(kk);
                
            else
                c=lamda_new(kk);
            end
        end
        if d-c<tol
            p_1_wf=p_1l_new(kk);
            p_2_wf=p_2l_new(kk);
            p_3_wf=p_3l_new(kk);
            break
        end
end
for n=1:length(m)
    p_11(n)=p_1_wf;
    p_21(n)=p_2_wf;
    p_31(n)=p_3_wf;
end

plot(m,p_1,'linewidth',1.25)
hold on
grid on
plot(m,p_2,':','linewidth',1.25)
plot(m,p_3,'-.','linewidth',1.25)

plot(m,p_11,'--','linewidth',1.25)
plot(m,p_21,'--','linewidth',1.25)
plot(m,p_31,'--','linewidth',1.25)

xlabel('log_2M')
ylabel('Power Allocation')
legend('g1=-20dB','g2=0dB','g3=20dB','WF')

text(6,53,'ch3')
text(6,47,'ch2')
text(6,2,'ch1')