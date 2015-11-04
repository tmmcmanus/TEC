function[h11_central,h12_central_x,hx_calc,f11_central] = partial_x(Z,X,i,j);
for n=2:(i.*j)-1
    h11_central(n)=1+(((Z(n+1)-Z(n-1)))./((X(n+1)-X(n-1)))).^2;
    h12_central_x(n)=(Z(n+1)-Z(n-1))./(X(n+1)-X(n-1));
    f11_central(n)=((Z(n+1)+Z(n-1)-2.*Z(n))./(X(n-1)-X(n)).^2);
end

%To be used with Z_virtual
h11_central(i:i:end)=[]; 
h11_central(1:i-1:end)=[];
h11_central(1:i-2)=[];
h11_central(end-(i-3):end)=[];

%To be used when using Z_virtual
 h12_central_x(i:i:end)=[];
 h12_central_x(1:i-1:end)=[];
 h12_central_x(1:i-2)=[];
 h12_central_x(end-(i-3):end)=[];
 hx_calc = h12_central_x;

%To be used with Z_virtual
 f11_central(i:i:end)=[];
 f11_central(1:i-1:end)=[];
 f11_central(1:i-2)=[];
 f11_central(end-(i-3):end)=[];
end
