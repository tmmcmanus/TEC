function[h11_central,h11_new,h12_central_x,f11_central]=new_partial_x(Z,X,i,j)
%INVESTIGATE - UNKNOWN IMPLEMENTATION 10-9-2013
A=norm_n_row3.';
grp=i;
numB=numel(A)+fix(numel(A)./grp);
idx=true(numB,1);
idx(grp:grp:numB)=false;
B(idx)=A;
B_new = zeros(length(Z),1);


for n = 1:length(B)
    B_new(n) = B(n);
end
Z_new = Z+B_new;

for n=2:(i.*j)-1
    h11_central(n)=1+(((Z_new(n+1)-Z_new(n-1)))./((X(n+1)-X(n-1)))).^2;
    h11_new(n)=1+(((Z_new(n+1)-Z_new(n-1)))./(2.*(i_virtual(n-1)+i_virtual(n)))).^2;
    h12_central_x(n)=(Z_new(n+1)-Z_new(n-1))./(X(n+1)-X(n-1));
    f11_central(n)=((Z_new(n+1)+Z_new(n-1)-2.*Z_new(n))./(X(n-1)-X(n)).^2);
end
%To be used in conjunction with Z_virtual_new
h11_central(i-1:i:end)=[];
h11_central(i:i-1:end)=[];
h11_central(1:i-2)=[];
h11_central(end-(i-3):end)=[];


%To be used when using Z_virtual_new
h12_central_x(i-1:i:end)=[];
h12_central_x(i:i-1:end)=[];
h12_central_x(1:i-2)=[];
h12_central_x(end-(i-3):end)=[];

%To be used with Z_virtual_new
 f11_central(i:i:end)=[];
 f11_central(1:i-1:end)=[];
 f11_central(1:i-2)=[];f11_central(end-(i-3):end)=[];
end
