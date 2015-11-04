function[pdx_gamma_2_12, pdx_gamma_1_21,pdx_gamma_2_22,pdx_gamma_1_22] = pdx_gamma(i,X_f2,gamma_2_12,gamma_1_21,gamma_2_22,gamma_1_22);
limit = length(gamma_2_12);
for n=2:limit-1;
    pdx_gamma_2_12(n) = (gamma_2_12(n+1)-gamma_2_12(n-1))./(X_f2(n+1)-X_f2(n-1));
    pdx_gamma_1_21(n) = (gamma_1_21(n+1)-gamma_1_21(n-1))./(X_f2(n+1)-X_f2(n-1));
    pdx_gamma_2_22(n) = (gamma_2_22(n+1)-gamma_2_22(n-1))./(X_f2(n+1)-X_f2(n-1));
    pdx_gamma_1_22(n) = (gamma_1_22(n+1)-gamma_1_22(n-1))./(X_f2(n+1)-X_f2(n-1));
end

%Stripping Stage
pdx_gamma_2_12(i-4:i-4:end)=[];
pdx_gamma_2_12(1:i-5:end)=[];
pdx_gamma_2_12(1:i-6)=[];
pdx_gamma_2_12(end-(i-7):end)=[];

pdx_gamma_1_21(i-4:i-4:end)=[];
pdx_gamma_1_21(1:i-5:end)=[];
pdx_gamma_1_21(1:i-6)=[];
pdx_gamma_1_21(end-(i-7):end)=[];

pdx_gamma_2_22(i-4:i-4:end)=[];
pdx_gamma_2_22(1:i-5:end)=[];
pdx_gamma_2_22(1:i-6)=[];
pdx_gamma_2_22(end-(i-7):end)=[];

pdx_gamma_1_22(i-4:i-4:end)=[];
pdx_gamma_1_22(1:i-5:end)=[];
pdx_gamma_1_22(1:i-6)=[];
pdx_gamma_1_22(end-(i-7):end)=[];

end
