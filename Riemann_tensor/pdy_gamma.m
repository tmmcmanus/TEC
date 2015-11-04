function[pdy_gamma_2_11,pdy_gamma_1_11,pdy_gamma_2_12,pdy_gamma_1_21]=pdy_gamma(i,Y_f2,gamma_2_11,gamma_1_11,gamma_2_12,gamma_1_21);
limit = length(gamma_2_11);
for n=1:limit-2*(i-4);
    pdy_gamma_2_11(n) = (gamma_2_11(2*(i-4)+n)-gamma_2_11(n))./(Y_f2(2*(i-4)+n)-Y_f2(n));
    pdy_gamma_1_11(n) = (gamma_1_11(2*(i-4)+n)-gamma_1_11(n))./(Y_f2(2*(i-4)+n)-Y_f2(n));
    pdy_gamma_2_12(n) = (gamma_2_12(2*(i-4)+n)-gamma_2_12(n))./(Y_f2(2*(i-4)+n)-Y_f2(n));
    pdy_gamma_1_21(n) = (gamma_1_21(2*(i-4)+n)-gamma_1_21(n))./(Y_f2(2*(i-4)+n)-Y_f2(n));
end


pdy_gamma_2_11(1:i-4:end)=[];
pdy_gamma_2_11(i-5:i-5:end)=[];

pdy_gamma_1_11(1:i-4:end)=[];
pdy_gamma_1_11(i-5:i-5:end)=[];

pdy_gamma_2_12(1:i-4:end)=[];
pdy_gamma_2_12(i-5:i-5:end)=[];

pdy_gamma_1_21(1:i-4:end)=[];
pdy_gamma_1_21(i-5:i-5:end)=[];
end