function[R_11,R_12,R_21,R_22,R] = Ricci_2nd_kind(R_2_121,R_1_112,R_1_212,contra_h11,contra_h12,contra_h22,i)
R_11 = R_2_121;
R_12 = R_1_112;
R_21 = R_12;
R_22 = R_1_212;

%Striping contraviariant tensors of order 2
contra_h11_strip = contra_h11';
contra_h11_strip(i-4:i-4:end)=[];
contra_h11_strip(1:i-5:end)=[];
contra_h11_strip(1:i-6)=[];
contra_h11_strip(end-(i-7):end)=[];

contra_h12_strip = contra_h12';
contra_h12_strip(i-4:i-4:end)=[];
contra_h12_strip(1:i-5:end)=[];
contra_h12_strip(1:i-6)=[];
contra_h12_strip(end-(i-7):end)=[];

contra_h21_strip = contra_h12_strip;

contra_h22_strip = contra_h22';
contra_h22_strip(i-4:i-4:end)=[];
contra_h22_strip(1:i-5:end)=[];
contra_h22_strip(1:i-6)=[];
contra_h22_strip(end-(i-7):end)=[];

R = contra_h11_strip.*R_11' + 2.*(contra_h12_strip).*(R_12') + contra_h22_strip.*R_22';
end