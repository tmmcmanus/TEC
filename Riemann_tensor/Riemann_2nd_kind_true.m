function[R_2_112_true,R_1_112_true,R_2_221_true,R_1_212_true,R_2_121_true,R_1_121_true,R_2_212_true,R_1_221_true,X_f3,Y_f3]=Riemann_2nd_kind_true(X_f2,Y_f2,gamma_1_21_true,gamma_2_11_true,gamma_2_12_true,gamma_1_22_true,i,j);

X_f3=X_f2;
X_f3(i-4:i-4:end)=[];
X_f3(1:i-5:end)=[];
X_f3(1:i-6)=[];
X_f3(end-(i-7):end)=[];

Y_f3=Y_f2;
Y_f3(i-4:i-4:end)=[];
Y_f3(1:i-5:end)=[];
Y_f3(1:i-6)=[];
Y_f3(end-(i-7):end)=[];

gamma_1_21_strip = gamma_1_21_true;
gamma_1_21_strip(i-4:i-4:end)=[];
gamma_1_21_strip(1:i-5:end)=[];
gamma_1_21_strip(1:i-6)=[];
gamma_1_21_strip(end-(i-7):end)=[];

gamma_2_11_strip = gamma_2_11_true;
gamma_2_11_strip(i-4:i-4:end)=[];
gamma_2_11_strip(1:i-5:end)=[];
gamma_2_11_strip(1:i-6)=[];
gamma_2_11_strip(end-(i-7):end)=[];

gamma_2_12_strip = gamma_2_12_true;
gamma_2_12_strip(i-4:i-4:end)=[];
gamma_2_12_strip(1:i-5:end)=[];
gamma_2_12_strip(1:i-6)=[];
gamma_2_12_strip(end-(i-7):end)=[];

gamma_2_21_strip = gamma_2_12_strip;

gamma_1_22_strip = gamma_1_22_true;
gamma_1_22_strip(i-4:i-4:end)=[];
gamma_1_22_strip(1:i-5:end)=[];
gamma_1_22_strip(1:i-6)=[];
gamma_1_22_strip(end-(i-7):end)=[];


alpha = (1-(X_f3).^2 + (Y_f3).^2)./((1+(X_f3).^2 + (Y_f3).^2).^2);
R_2_112_true = alpha + ((gamma_1_21_strip).*(gamma_2_11_strip)) + ((gamma_2_12_strip).*(gamma_2_21_strip));
R_2_121_true = -1.*R_2_112_true;

alpha = (-2.*(X_f3).*(Y_f3))./((1+(X_f3).^2 + (Y_f3).^2).^2);
R_1_112_true = alpha + (gamma_2_12_strip.*gamma_1_21_strip);
R_1_121_true = -1.*R_1_112_true;

alpha = (-2.*(X_f3).*(Y_f3))./((1+(X_f3).^2 + (Y_f3).^2).^2);
R_2_221_true = alpha + (gamma_1_21_strip).*(gamma_2_12_strip);
R_2_212_true = -1.*R_2_221_true;

alpha = (-1-(X_f3).^2 + (Y_f3).^2)./((1+(X_f3).^2 + (Y_f3).^2).^2);
R_1_212_true = alpha - (gamma_1_21_strip).*(gamma_1_21_strip) - (gamma_2_21_strip).*(gamma_1_22_strip);
R_1_221_true = -1.*R_1_212_true;
end