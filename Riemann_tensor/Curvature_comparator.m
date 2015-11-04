function[K,K_R_error_true,K_R_error] = Curvature_comparator(R,R_true,X_f3,Y_f3)
%Assuming map comprised of Monge Patches z = h(u,v)
K =-1./((1+(X_f3).^2 + (Y_f3).^2)).^2;
%Recalling that R = 2*K;
K_R_error_true = abs(K - (R_true./2));
K_R_error = abs(K - (R./2));
end