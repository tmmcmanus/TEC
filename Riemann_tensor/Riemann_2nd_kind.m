function[R_2_112,R_1_112,R_2_221,R_1_212,R_2_121,R_1_121,R_2_212,R_1_221] = Riemann_2nd_kind(gamma_2_12,gamma_1_21,gamma_2_22,gamma_1_22,gamma_2_11,gamma_1_11,i,j,pdx_gamma_2_21,pdy_gamma_2_11,pdx_gamma_1_21,pdy_gamma_1_11,pdy_gamma_2_12,pdx_gamma_2_22,pdx_gamma_1_22,pdy_gamma_1_21);
gamma_2_12(1:i-4)=[];
gamma_2_12(i-4:i-4:end)=[];
gamma_2_12(1:i-5:end)=[];
gamma_2_12(end-(i-7):end)=[];

gamma_1_21(1:i-4)=[];
gamma_1_21(i-4:i-4:end)=[];
gamma_1_21(1:i-5:end)=[];
gamma_1_21(end-(i-7):end)=[];

gamma_2_22(1:i-4)=[];
gamma_2_22(i-4:i-4:end)=[];
gamma_2_22(1:i-5:end)=[];
gamma_2_22(end-(i-7):end)=[];

gamma_1_22(1:i-4)=[];
gamma_1_22(i-4:i-4:end)=[];
gamma_1_22(1:i-5:end)=[];
gamma_1_22(end-(i-7):end)=[];

gamma_2_11(1:i-4)=[];
gamma_2_11(i-4:i-4:end)=[];
gamma_2_11(1:i-5:end)=[];
gamma_2_11(end-(i-7):end)=[];

gamma_1_11(1:i-4)=[];
gamma_1_11(i-4:i-4:end)=[];
gamma_1_11(1:i-5:end)=[];
gamma_1_11(end-(i-7):end)=[];


R_1_111 = zeros((i-6).*(j-6),1);
R_1_211 = R_1_111;
R_1_122 = R_1_111;
R_1_222 = R_1_111;
R_2_111 = zeros((i-6).*(j-6),1);
R_2_211 = R_2_111;
R_2_122 = R_2_111;
R_2_222 = R_2_111;

R_2_112 = pdx_gamma_2_21 - pdy_gamma_2_11 + (gamma_1_21).*(gamma_2_11) - (gamma_1_11).*(gamma_2_12) + (gamma_2_12).*(gamma_2_12) - (gamma_2_11).*(gamma_2_22);
R_1_112 = pdx_gamma_1_21 - pdy_gamma_1_11 + (gamma_2_12).*(gamma_1_21) - (gamma_2_11).*(gamma_1_22);
R_2_221 = pdy_gamma_2_12 - pdx_gamma_2_22 + (gamma_1_21).*(gamma_2_12) - (gamma_1_22).*(gamma_2_11);
R_1_212 = pdx_gamma_1_22 - pdy_gamma_1_21 + (gamma_1_22).*(gamma_1_11) - (gamma_1_21).*(gamma_1_21) + (gamma_2_22).*(gamma_1_21) - (gamma_2_12).*(gamma_1_22);

R_2_121 = -1*R_2_112;
R_1_121 = -1*R_1_112;
R_2_212 = -1*R_2_221;
R_1_221 = -1*R_1_212;
end
