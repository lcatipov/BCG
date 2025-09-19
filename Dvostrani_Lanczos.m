function [V,W,T]=Dvostrani_Lanczos(A,maxiter)
  
  n=size(A,1);
  v1=zeros(n,1);      %dani vektor v1 t.d. norm(v1,2)=1
  w1=zeros(n,1);        %dani vektor w1 t.d. <v1,w1>=1;
  v1(1)=1;
  w1(1)=1;
  
  %v1=v1/norm(v1,2);
  %w1=w1/norm(v1'*w1);
  
  beta=0; gama=0;
  Alpha=[];
  Beta=[];
  Gama=[];
  V=[];
  W=[];
  
  v0=zeros(n,1);
  w0=zeros(n,1);
  
  for k=1:maxiter
    V=[V,v1];
    W=[W,w1];
    Av=A*v1;
    Aw=A'*w1;
    alpha=Av'*w1;
    Alpha=[Alpha;alpha];
    vtilda=Av-alpha*v1-beta*v0;
    wtilda=Aw-conj(alpha)*w1-gama*w0;
    gama=norm(vtilda,2);
    
    if gama==0
      break;
    endif
    Gama=[Gama;gama];
    
    v1=vtilda/gama;
    beta=v1'*wtilda;
    
    if beta==0
      break;
    endif
    Beta=[Beta;beta];
    
    w1=wtilda/conj(beta);
    v0=V(:,end);
    w0=W(:,end);
    
  endfor
  vtilda
  T=diag(Alpha)+diag(Beta(1:end-1),1)+diag(Gama(1:end-1),-1);
end