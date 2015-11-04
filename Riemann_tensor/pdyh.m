function [pdyh11,pdyh12,pdyh22,Y_f1,pdyh11_raw,pdyh12_raw,pdyh22_raw] = pdyh(Y,j,i,h11_central,h12_central,h22_central);
Y_f1 = Y;
Y_f1(1:j)=[]; 
Y_f1(1:j:end)=[]; 
Y_f1(j-1:j-1:end)=[]; 
Y_f1(end-(j-3):end)=[];

limit = length(h11_central);

for n=1:limit - 2*(i-2);
    pdyh11(n) = (h11_central(2*(i-2)+n)-h11_central(n))./(Y_f1(2*(i-2)+n)-Y_f1(n));
    pdyh12(n) = (h12_central(2*(i-2)+n)-h12_central(n))./(Y_f1(2*(i-2)+n)-Y_f1(n));
    pdyh22(n) = (h22_central(2*(i-2)+n)-h22_central(n))./(Y_f1(2*(i-2)+n)-Y_f1(n));
end

pdyh11_raw = pdyh11;
pdyh12_raw = pdyh12;
pdyh22_raw = pdyh22;

pdyh11(1:i-2:end)=[];
pdyh11(i-3:i-3:end)=[];

pdyh12(1:i-2:end)=[];
pdyh12(i-3:i-3:end)=[];

pdyh22(1:i-2:end)=[];
pdyh22(i-3:i-3:end)=[];
end