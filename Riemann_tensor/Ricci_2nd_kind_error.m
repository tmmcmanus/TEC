function[R_true,R_error] = Ricci_2nd_kind_error(R_2_121_true, R_1_112_true,R_1_212_true,X_f3,Y_f3,R)
R_11_true = R_2_121_true;
R_12_true = R_1_112_true;
R_21_true = R_12_true;
R_22_true = R_1_212_true;

alpha = 1./((((1+(Y_f3).^2).*(1+(X_f3).^2)) - ((X_f3).^2).*(Y_f3).^2));
R_true =  alpha.*(1+(X_f3.^2)).*R_11_true + 2.*(alpha).*(-1.*(X_f3).*(Y_f3)).*R_12_true + alpha.*(1+(Y_f3).^2).*R_22_true;
R_error = abs(R_true - R);
end