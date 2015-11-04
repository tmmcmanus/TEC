function[f12_central] = xy_partial(Z,X,Y,i,j)
for n=i+2:(i-1).*j-1;
    nume(n)=(Z(n+(i+1))-Z(n+(i-1))-Z(n-(i-1))+Z(n-(i+1)));
    denom(n)=(4.*(((X(n+1)-X(n)).*(Y(n+i)-Y(n)))));
    f12_central(n)=nume(n)./denom(n);
end

%removing first layer
f12_central(1:i+1)=[];
% Z_values(1:i+1)=[];
% denom(1:i+1)=[];
% nume(1:i+1)=[];
%removing discontinities along the u
f12_central(i:i:end)=[];
% Z_values(i:i:end)=[];
% denom(i:i:end)=[];
% nume(i:i:end)=[];
f12_central(i-1:i-1:end)=[];
% Z_values(i-1:i-1:end)=[];
% denom(i-1:i-1:end)=[];
% nume(i-1:i-1:end)=[];
end