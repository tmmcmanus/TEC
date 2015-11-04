function[gamma_1_11,gamma_1_21,gamma_1_22,gamma_2_11,gamma_2_12,gamma_2_22,gamma_1_11_true,gamma_1_21_true,gamma_1_22_true,gamma_2_11_true,gamma_2_12_true,gamma_2_22_true,contra_h11,contra_h12,contra_h22] = Chris_second_kind(i,j,h11_central,h12_central,h22_central,gamma_111,gamma_112,gamma_221,gamma_222,gamma_121,gamma_122,X_f2,Y_f2,gamma_111_true,gamma_121_true,gamma_221_true,gamma_112_true,gamma_122_true,gamma_222_true);

contra = cell(1,length(h11_central));
for n=1:length(h11_central)
    contra{n}=inv([h11_central(n),h12_central(n);h12_central(n),h22_central(n)]);
    contra_h11(n) = contra{n}(1);
    contra_h12(n) = contra{n}(2);
    contra_h22(n) = contra{n}(4);
end

%Stripping procedure
contra_h11(1:i-2)=[];
contra_h11(i-2:i-2:end)=[];
contra_h11(1:i-3:end)=[];
contra_h11(end-(i-5):end)=[];

contra_h12(1:i-2)=[];
contra_h12(i-2:i-2:end)=[];
contra_h12(1:i-3:end)=[];
contra_h12(end-(i-5):end)=[];

contra_h22(1:i-2)=[];
contra_h22(i-2:i-2:end)=[];
contra_h22(1:i-3:end)=[];
contra_h22(end-(i-5):end)=[];

%Populating known analytic solution
coeff = 1./((1+(Y_f2).^2).*(1+(X_f2).^2) - ((X_f2.^2).*(Y_f2.^2)));
contra_h11_true = (1+X_f2.^2).*coeff;
contra_h12_true = -1.*((X_f2).*(Y_f2)).*coeff;
contra_h22_true = (1+Y_f2.^2).*coeff;
gamma_1_11 = contra_h11.*gamma_111 + contra_h12.*gamma_112;
gamma_1_21 = contra_h11.*gamma_121 + contra_h12.*gamma_122;
gamma_1_22 = contra_h11.*gamma_221 + contra_h12.*gamma_222;
gamma_2_11 = contra_h12.*gamma_111 + contra_h22.*gamma_112;
gamma_2_12 = contra_h12.*gamma_121 + contra_h22.*gamma_122;
gamma_2_22 = contra_h12.*gamma_221 + contra_h22.*gamma_222;


gamma_1_11_true = contra_h11_true.*gamma_111_true + contra_h12_true.*gamma_112_true;
gamma_1_21_true = contra_h11_true.*gamma_121_true + contra_h12_true.*gamma_122_true;
gamma_1_22_true = contra_h11_true.*gamma_221_true + contra_h12_true.*gamma_222_true;
gamma_2_11_true = contra_h12_true.*gamma_111_true + contra_h22_true.*gamma_112_true;
gamma_2_12_true = contra_h12_true.*gamma_121_true + contra_h22_true.*gamma_122_true;
gamma_2_22_true = contra_h12_true.*gamma_221_true + contra_h22_true.*gamma_222_true;
end