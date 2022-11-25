function M = EnsambleMasa1D(x,r)
% Matriz de masa para problema de Dirichlet homogeneo.

n = length(x)-2;
M = zeros(n,n);  
M(1,1)=  r( 0.5*(x(2)+ x(1)) )* ( x(2)- x(1))/3;
M(n,n)=  r( 0.5*(x(n+2)+ x(n+1)) )*  ( x(n+2) - x(n+1))/3;
for i = 2:n 
    h = x(i+1) - x(i); % longitud del intervalo
    xmid = (x(i+1) + x(i))/2; % punto medio del intervalo
    rmid = r(xmid); % valor de r(x) en el punto medio

    M(i-1,i-1) = M(i-1,i-1) +  rmid*h/3; % suma rmid*h/3 a M(i,i)
    M(i-1,i) = M(i-1,i) +  rmid*h/6;
    M(i,i-1) = M(i,i-1) + rmid*h/6;
    M(i,i) = M(i,i) + rmid*h/3;
end

