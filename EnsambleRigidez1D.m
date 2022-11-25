function A = EnsambleRigidez1D(x,p)
% Matriz de rigidez para problema de Dirichlet homogeneo.
  
n = length(x)-2;   % x0< x1 < ... < xN n+1 puntos en la malla %
A = zeros(n,n);
A(1,1)= p( 0.5*(x(2)+ x(1)) )/ ( x(2)- x(1));
A(n,n)= p( 0.5*(x(n+2)+ x(n+1)) )/ ( x(n+2) - x(n+1));
for i = 2:n
    h = x(i+1) - x(i);
    xmid = (x(i+1) + x(i))/2; % punto medio del intervalo
    pmid = p(xmid); % valor de p(x) en el punto medio
    A(i-1,i-1) = A(i-1,i-1) + pmid/h; % suma pmid/h a A(i-1,i-1)
    A(i-1,i) = A(i-1,i) - pmid/h;
    A(i,i-1) = A(i,i-1) - pmid/h;
    A(i,i) = A(i,i) + pmid/h;
end

