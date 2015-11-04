function [pdxh11,pdxh12,pdxh22,X_f1] = pdxh(X,i,h11_central,h12_central,h22_central);
X_f1 = X;
X_f1(1:i)=[]; 
X_f1(1:i:end)=[]; 
X_f1(i-1:i-1:end)=[]; 
X_f1(end-(i-3):end)=[];

limit = length(h11_central);
for n=2:limit-1
    pdxh11(n)=(h11_central(n+1)-h11_central(n-1))./(X_f1(n+1)-X_f1(n-1));
    pdxh12(n)=(h12_central(n+1)-h12_central(n-1))./(X_f1(n+1)-X_f1(n-1));
    pdxh22(n)=(h22_central(n+1)-h22_central(n-1))./(X_f1(n+1)-X_f1(n-1));
end
pdxh11(i-2:i-2:end)=[];
pdxh11(1:i-3:end)=[];
pdxh11(1:i-4)=[];
pdxh11(end-(i-5):end)=[];

pdxh12(i-2:i-2:end)=[];
pdxh12(1:i-3:end)=[];
pdxh12(1:i-4)=[];
pdxh12(end-(i-5):end)=[];

pdxh22(i-2:i-2:end)=[];
pdxh22(1:i-3:end)=[];
pdxh22(1:i-4)=[];
pdxh22(end-(i-5):end)=[];

end
