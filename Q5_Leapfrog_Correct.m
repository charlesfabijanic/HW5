clc
clearvars

h_x=0.001;
h_t=0.0001;
x02=0:h_x:1;
% x1=0.2+h_x:h_x:1;

length_x=length(x02);%+length(x1);
t=linspace(h_t,8,8/h_t);
length_t=length(t);
u=0.08;

T=zeros(length_t,length_x);
%introducing BC


%introducing IC
x=x02;
for i=1:length_x
    if x(i)<=0.2
        T(1,i)=1-(10*x(i)-1)^2;
    else
        T(1,i)=0;
    end
end
T(:,1)=0;
T(:,end)=0;
%solving matrix Explicit Euler time advancement and second-order central
%differencing

clearvars i
g=u*h_t/(2*h_x);
for j=2:length_x-1
    T(2,j)=T(1,j)-g*(T(1,j+1)-T(1,j-1));
end
g=u*h_t/(h_x);
for i=2:length_t-1
    for j=2:length_x-1
        T(i+1,j)=T(i-1,j)-g*(T(i,j+1)-T(i,j-1));
    end
end

xx=0:0.005:1;
tt=0:4:8;
for i=1:length(tt)
    for j=1:length(xx)
        if (xx(j)-u*tt(i))<=0.2
            TT(i,j)=1-(10*(xx(j)-u*tt(i))-1)^2;
        else
            TT(i,j)=0;
        end
    end
end
ind=find(TT<=0);
TT(ind)=0;
ind2=linspace(1,length_t,3);
ind2=round(ind2);
plot(x,T(ind2,:))
hold on
plot(xx,TT,'--')
legend('LF t=0','LF t=4','LF t=8','Exact t=0','Exact t=4','Exact t=8')