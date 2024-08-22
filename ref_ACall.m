function [res,err]=ref_ACall(opt,S0,k,r,sigma,T,M,npt)

if(opt==1)
   dt=T/M;
   a=(r-sigma^2/2)*dt;
   b=sigma*sqrt(dt);

   ff1=0;
   ff2=0;

   for i=1:npt
       stm=zeros(M,1);
       xi=ournormal;
       stm(1)=S0*exp(a+b*xi(1));
       for j=2:M
           xi=ournormal;
           stm(j)=stm(j-1)*exp(a+b*xi(1));
       end
       tmp=exp(-r*T)*max([(sum(stm)/M-k) 0]);
       ff1=ff1+tmp;
       ff2=ff2+tmp^2;
   end

   res=ff1/npt;
   err=sqrt((ff2/npt-res^2)/npt);


elseif(opt==2)
   dt=T/M;
   a=(r-sigma^2/2)*dt;
   b=sigma*sqrt(dt);

% this only for comparison
   nn=round(npt/2);

   ff1=0;
   ff2=0;

   for i=1:nn
       stmp=zeros(M,1);
       stmm=zeros(M,1);  
       xi=ournormal;
       stmp(1)=S0*exp(a+b*xi(1));
       stmm(1)=S0*exp(a-b*xi(1));
       for j=2:M
           xi=ournormal;
           stmp(j)=stmp(j-1)*exp(a+b*xi(1));
           stmm(j)=stmm(j-1)*exp(a-b*xi(1));
       end

       tmp=exp(-r*T)*(max([(sum(stmp)/M-k) 0])+max([(sum(stmm)/M-k) 0]))/2;
       ff1=ff1+tmp;
       ff2=ff2+tmp^2;
   end

   res=ff1/nn;
   err=sqrt((ff2/nn-res^2)/nn);



elseif(opt==3)

   dt=T/M;
   a=(r-sigma^2/2)*dt;
   b=sigma*sqrt(dt);

% first find covariance and variate variance
   nvars=1e4;
   phis=zeros(nvars,1);
   cvs =zeros(nvars,1);

   for i=1:nvars
       stm=zeros(M,1);
       xi=ournormal;
       stm(1)=S0*exp(a+b*xi(1));
       for j=2:M
           xi=ournormal;
           stm(j)=stm(j-1)*exp(a+b*xi(1));
       end
       cvs(i) =exp(-r*M*dt)*(sum(stm)+S0)/M;
       phis(i)=exp(-r*T)*max([(sum(stm)/M-k) 0]);
   end
   mucvs =sum(cvs )/nvars;
   muphis=sum(phis)/nvars;
   mu2cvs =sum(cvs.*cvs )/nvars;   

   cvsvar=mu2cvs-mucvs^2;

   covar=sum((cvs-mucvs).*(phis-muphis))/nvars;
   bval=covar/cvsvar;
% this only for comparison
   nn=npt;

   meancv=S0*exp(-r*dt*M)*(exp((M+1)*r*dt) -1)/(M*(exp(r*dt) -1));

   
   ff1=0;
   ff2=0;

   for i=1:nn
       stm=zeros(M,1);
       xi=ournormal;
       stm(1)=S0*exp(a+b*xi(1));
       for j=2:M
           xi=ournormal;
           stm(j)=stm(j-1)*exp(a+b*xi(1));
       end
       cv =exp(-r*M*dt)*(sum(stm)+S0)/M;
       tmp=exp(-r*T)*max([(sum(stm)/M-k) 0])-bval*(cv-meancv);
       ff1=ff1+tmp;
       ff2=ff2+tmp^2;
   end

   res=ff1/nn;
   err=sqrt((ff2/nn-res^2)/nn);

elseif(opt==4)

   dt=T/M;
   a=(r-sigma^2/2)*dt;
   b=sigma*sqrt(dt);

% first find covariance and variate variance
   nvars=1e4;
   phis=zeros(nvars,1);
   cvs =zeros(nvars,1);

   for i=1:nvars
       stm=zeros(M,1);
       xi=ournormal;
       stm(1)=S0*exp(a+b*xi(1));
       for j=2:M
           xi=ournormal;
           stm(j)=stm(j-1)*exp(a+b*xi(1));
       end
       cvs(i) =exp(-r*T)*max([(prod(stm)^(1/M)-k) 0]);
       phis(i)=exp(-r*T)*max([(sum(stm)/M-k) 0]);
   end
   mucvs =sum(cvs )/nvars;
   muphis=sum(phis)/nvars;
   mu2cvs =sum(cvs.*cvs )/nvars;   

   cvsvar=mu2cvs-mucvs^2;

   covar=sum((cvs-mucvs).*(phis-muphis))/nvars;
   bval=covar/cvsvar;
% this only for comparison
   nn=npt;

   meancv=GeometricAsian (S0, k ,r ,T, sigma,0,M);

%   error
   
   ff1=0;
   ff2=0;

   for i=1:nn
       stm=zeros(M,1);
       xi=ournormal;
       stm(1)=S0*exp(a+b*xi(1));
       for j=2:M
           xi=ournormal;
           stm(j)=stm(j-1)*exp(a+b*xi(1));
       end
       cv =exp(-r*T)*max([(prod(stm)^(1/M)-k) 0]);
       tmp=exp(-r*T)*max([(sum(stm)/M-k) 0])-bval*(cv-meancv);
       ff1=ff1+tmp;
       ff2=ff2+tmp^2;
   end

   res=ff1/nn;
   err=sqrt((ff2/nn-res^2)/nn);

elseif(opt==5)
   dt=T/M;
   a=(r-sigma^2/2)*dt;

   ns=20; % number of strata
   npts=round(npt/ns);

   res=0;
   err=0;
   for j=1:ns
      ff1=0;
      ff2=0;
      for i=1:npts
        stm=zeros(M,1);
        z=norminv(rand(1)/ns+(j-1)/ns);
        wtm=z*sqrt(T); % wiener increment at the last step
        ww=BrownianBridge(M,dt,wtm);
        % stratify the last increment
        stm(1)=S0*exp(a+sigma*ww(1)); % w already contains dt factor
        for j1=2:M
            stm(j1)=stm(j1-1)*exp(a+sigma*ww(j1));
        end
        tmp=exp(-r*T)*max([(sum(stm)/M-k) 0])/ns;
        ff1=ff1+tmp;
        ff2=ff2+tmp^2;
     end

    res=res+ff1/npts;
    err=err+(ff2/npts-(ff1/npts)^2)/npts;
   end
   err=sqrt(err);

end 

end