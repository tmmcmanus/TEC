function[L_det,M_det,N_det,K_det]=algebraic_SFF(f11_central,h11_central,h22_central,h12_central,f12_central,f22_central)
L_det=zeros(length(f11_central),1);
M_det=zeros(length(f11_central),1);
N_det=zeros(length(f11_central),1);
for n=1:length(f11_central)
    L_det(n)=(f11_central(n))./(sqrt(h11_central(n).*h22_central(n) - (h12_central(n)).^2));
    M_det(n)=(f12_central(n))./(sqrt(h11_central(n).*h22_central(n) - (h12_central(n)).^2));
    N_det(n)=(f22_central(n))./(sqrt(h11_central(n).*h22_central(n) - (h12_central(n)).^2));
    K_det(n)= (L_det(n).*N_det(n) - (M_det(n)).^2)./(h11_central(n).*h22_central(n) - (h12_central(n)).^2);
end
end