function [x,k,rez,rezh]=BCG(A,b,tol,maxiter)
  
  n=size(A,1);
  x=zeros(n,1);
  
  r=b-A*x;
  rh=r;
  
  if r'*rh==0
    error;
  endif
  
  p=r;
  ph=rh;
  
  for k=1:maxiter
    Ap=A*p;
    Aph=A'*ph;
    alpha=(r'*rh)/(Ap'*ph);
    x=x+alpha*p;
    r1=r-alpha*Ap;
    if norm(r1)<tol
      break;
    endif
    r1h=rh-conj(alpha)*Aph;
    if norm(r1h)<tol
      break;
    endif
    beta=(r1'*r1h)/(r'*rh);
    p1=r1+beta*p;
    p1h=r1h+conj(beta)*ph;
    
    r=r1;
    rh=r1h;
    p=p1;
    ph=p1h;
  endfor
  
  rez=norm(r);
  rezh=norm(rh);
  
end