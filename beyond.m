clear all
close all
clc
a=-1; b=1;                      % spatial interval
N=1000; h=(b-a)/N;
x=a:h:b;
h1=(b-a)/(N-2);
y=a:h1:b;
p=@(x) 2;                       % p:= sigma^2/(2R^2)
r=@(x) 1;
ni=12; dt=0.1;                  % temporal interval

A=EnsambleRigidez1D(x,p);       % Assembled matrices
M=EnsambleMasa1D(x,r);
s=zeros(length(x)-2,ni);

for i=1:length(x)-2             % Initial condition
    if x(i)<0
        s(i,1)=0;
    else
        s(i,1)=1;
    end
end

tic
figure()                        % Solving the equation
for i=2:ni
    s(:,i) = ( M + dt*A )\(M*s(:,i-1));%+M*(s(:,i-1).*(1-s(:,i-1))));
    %show(x,full(s(:,i)));
    %s(:,i).*(1-(s:,i));
    plot(y,s(:,i))
    xlabel('z (number associated with each person)')
    ylabel('productivity')
    pause(0.1)
end




