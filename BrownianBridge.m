function dw=BrownianBridge(M,dt,wtm)
dw=zeros(M,1);
ww=zeros(M,1);
xi=ournormal;
mui=wtm/M; % dt cancels out
sigi=sqrt((M-1)*dt/M);
ww(1)=mui+sigi*xi(1);
dw(1)=ww(1);
for i=2:M-1
    xi=ournormal;
    mui=((M-i)*ww(i-1)+wtm)/(M-i+1); % dt cancels out
    sigi=sqrt((M-i)*dt/(M-i+1));
    ww(i)=mui+sigi*xi(1);
    dw(i)=ww(i)-ww(i-1);
end
dw(M)=wtm-ww(M-1);
end

 
