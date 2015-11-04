function[h22_central,h12_central_y,f22_central]=partial_y(Z,Y,i,j)

h22_central=zeros(1,(i.*(j-1)));
for n=i+1:(i.*(j-1))
    h22_central(n)=1+(((Z(n+i)-Z(n-i)))./((Y(n+i)-Y(n-i)))).^2;
    %Note change in h12_central_y 
    h12_central_y(n)=(Z(n+i)-Z(n-i))./(Y(n+i)-Y(n-i));
    %h12_central_y_2(n)=(Z(n+i)-Z(n-i))./(2.*(Y(n+i)-Y(n)));
    f22_central(n)=(Z(n+i)+Z(n-i)-2.*Z(n))./(Y(n+i)-Y(n)).^2;
end

h22_central(1:i)=[];
h22_central(1:i:end)=[];
h22_central(i-1:i-1:end)=[];

f22_central(1:i)=[];
f22_central(1:i:end)=[];
f22_central(i-1:i-1:end)=[];

h12_central_y(1:i)=[];
h12_central_y(1:i:end)=[];
h12_central_y(i-1:i-1:end)=[];

end