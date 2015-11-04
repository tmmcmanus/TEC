%TEC_2D
%Timothy M. McManus Jr.

clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing Pointwise .wrl file      %
% and extracting desired cooridnates %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load Virtual Space
fid = fopen('Shoe_Surface_129x129.wrl','r'); % Pointwise export file
textscan(fid, '%*s %*[^\n]', 5); % Skip header text
Virtual_data = textscan(fid, '%f %f %f %*[^\n]', 'delimiter','','MultipleDelimsAsOne', 1, 'CollectOutput', 1); % read numeric data
fclose(fid);
X_virtual = Virtual_data{1}(:,1);
Y_virtual = Virtual_data{1}(:,2);
Z_virtual = Virtual_data{1}(:,3);
Virtual=[X_virtual,Y_virtual,Z_virtual];

% %Load Virtual Permittivity Map
% fid = fopen('','r'); %Handmade Perm map
% Virtual_perm_map = textscan(fid, '%f %f %f %f'); %read virtual perm map
% fclose(fid);
% virtual_eps_11 = Virtual_perm_map{1};
% virtual_eps_12 = Virtual_perm_map{2};
% virtual_eps_21 = Virtual_perm_map{3};
% virtual_eps_22 = Virtual_perm_map{4};
% Virtual_perm = [virtual_eps_11,virtual_eps_12,virtual_eps_21,virtual_eps_22];

%Load Physical Space
fid = fopen('Flat_Plane_129x129.wrl','r');
textscan(fid,'%*s %*[^\n]', 5);
Physical_data = textscan(fid, '%f %f %f %*[^\n]', 'delimiter','','MultipleDelimsAsOne', 1, 'CollectOutput', 1);
fclose(fid);
X_physical = Physical_data{1}(:,1);
Y_physical = Physical_data{1}(:,2);
Z_physical = Physical_data{1}(:,3);
Physical=[X_physical,Y_physical,Z_physical];
%ensure correct orientation exist for both spaces

%Ensuring that the X_virtual_spacing and the X_physical_spacing are the
%same

%X_physical = X_virtual;

i=129; %Points in the i-direction 
j=129; %number of POINTS in the j-direction
eps_virtual=1;

%determining the number of vectors given in the i,j and k directions
i_vectors = (i-1).*(j);
j_vectors = (i).*(j-1);
cells=(i-1).*(j-1);

%Loading the Points that will be inputted into the FFF:
% X_fff=X_physical; Y_fff=Y_physical; Z_fff=Z_physical;
% 
% %Removing the Edges from X_fff,Y_fff and Z_fff s.t. the metric tensors can
% %be compared with those determined by the following method
% X_fff(i:i:end)=[]; X_fff(((i-1).*(j-1))+1:((i-1).*j))=[];
% Y_fff(i:i:end)=[]; Y_fff(((i-1).*(j-1))+1:((i-1).*j))=[];
% Z_fff(i:i:end)=[]; Z_fff(((i-1).*(j-1))+1:((i-1).*j))=[];
% 
% %Inputting X_fff, Y_fff and Z_fff into the Correct Parametric First
% %Fundamental Forms
% 
% %Hyperbolic Paraboloid
% E = 1+(Y_fff).^2; F=(X_fff.*Y_fff); G = 1+(X_fff).^2;
% 
% %%Prolate Spheroid
% % a = 1; b = 2; c = 3; v = acos(Z_fff./c); u = atan((a.*Y_fff)./(b.*X_fff));
% % E = (b.^2.*((cos(u)).^2) + a.^2.*((sin(u)).^2)).*(sin(u)).^2;
% % F = (b.^2 - a.^2).*(cos(u).*sin(u).*cos(v).*sin(v));
% % G = (a.^2.*(cos(u)).^2 +b.^2.*(sin(u)).^2).*(cos(v)).^2 + c.^2.*(sin(v)).^2;
% 
% %%Shoe Surface
%  %E = 1 + (X_fff).^4; F = -1.*((X_fff).^2).*(Y_fff); G = 1 + (Y_fff).^2;
% 
% %%Monkey Saddle
% %E = 1 +9.*(((X_fff).^2 - (Y_fff).^2)).^2; F = -18.*(X_fff.*Y_fff).*((X_fff).^2 - (Y_fff).^2); G = 1 + 36.*((X_fff.^2).*(Y_fff).^2);
% 
% for n=1:length(E)
%     Parametric_FFF_coefficients(n,:)=[E(n),F(n),G(n)];
% end
% 
% 
% %Populating the h tensor {New Approach}
% for n = 1:(i.*j)-1
%     h11(n) = 1 + ((Z_physical(n)-Z_physical(n+1))./(X_physical(n)-X_physical(n+1))).^2;
%     h12_x(n) = (Z_physical(n)-Z_physical(n+1))./(X_physical(n)-X_physical(n+1));
% end
% h11(i:i:end)=[]; h11(((i-1).*(j-1))+1:((i-1).*j))=[];
% h12_x(i:i:end)=[]; h12_x(((i-1).*(j-1))+1:((i-1).*j))=[];
% for n = 1:(i.*(j-1))
%     h22(n) = 1 + ((Z_physical(n)-Z_physical(n+i))./(Y_physical(n)-Y_physical(n+i))).^2;
%     h12_y(n)= (Z_physical(n)-Z_physical(n+i))./(Y_physical(n)-Y_physical(n+i));
% end
% % h22(i:i:end)=[];
% % h12_y(i:i:end)=[];
% % 
% % h12 = h12_x.*h12_y;
% % 
% % for n = 1:(i.*j)-1
% %     h11_new(n) = 1 + ((Y_physical(n)-Y_physical(n+1))./(X_physical(n)-X_physical(n+1))).^2 + ((Z_physical(n)-Z_physical(n+1))./(X_physical(n)-X_physical(n+1))).^2;
% %     h12_C_x_new(n) = (Z_physical(n)-Z_physical(n+1))./(X_physical(n)-X_physical(n+1));
% %     h12_B_x_new(n) = (Y_physical(n)-Y_physical(n+1))./(X_physical(n)-X_physical(n+1));
% % end
% % h11_new(i:i:end)=[]; h11_new(((i-1).*(j-1))+1:((i-1).*j))=[];
% % h12_C_x_new(i:i:end)=[]; h12_C_x_new(((i-1).*(j-1))+1:((i-1).*j))=[];
% % h12_B_x_new(i:i:end)=[]; h12_B_x_new(((i-1).*(j-1))+1:((i-1).*j))=[];
% % for n = 1:(i.*(j-1))
% %     h22_new(n) = ((X_physical(n)-X_physical(n+i))./(Y_physical(n)-Y_physical(n+i))).^2 + 1 + ((Z_physical(n)-Z_physical(n+i))./(Y_physical(n)-Y_physical(n+i))).^2;
% %     h12_C_y_new(n)= (Z_physical(n)-Z_physical(n+i))./(Y_physical(n)-Y_physical(n+i));
% %     h12_A_y_new(n) = (X_physical(n)-X_physical(n+i))./(Y_physical(n)-Y_physical(n+i));
% % end
% % h22_new(i:i:end)=[];
% % h12_C_y_new(i:i:end)=[];
% % h12_A_y_new(i:i:end)=[];
% % 
% % h12_C = h12_C_x_new.*h12_C_y_new;
% % h12_A = h12_A_y_new; h12_B = h12_B_x_new;
% % h12_new = h12_A + h12_B + h12_C;
% % 
% % %Percent Error Calculation
% % for n=1:length(h11_new)
% %     h11_error(n)= (abs(h11_new(n)-E(n))./abs(E(n))).*100;
% %     h22_error(n)= (abs(h22_new(n)-G(n))./abs(G(n))).*100;
% %     %note change to h12_error
% %     h12_error(n)= (abs(h12(n)-F(n))./abs(F(n))).*100;
% % end

for n = 1:(i.*j)-1
    a_physical(n)=sqrt((X_physical(n)-X_physical(n+1)).^2 + (Y_physical(n)-Y_physical(n+1)).^2 + (Z_physical(n) - Z_physical(n+1)).^2);
    a_virtual(n)= sqrt((X_virtual(n)-X_virtual(n+1)).^2 + (Y_virtual(n)-Y_virtual(n+1)).^2 + (Z_virtual(n) - Z_virtual(n+1)).^2);
    a_virtual_vector(n,:)=[(X_virtual(n)-X_virtual(n+1)),(Y_virtual(n)-Y_virtual(n+1)),(Z_virtual(n) - Z_virtual(n+1))];
    a_physical_vector(n,:)=[(X_physical(n)-X_physical(n+1)),(Y_physical(n)-Y_physical(n+1)),(Z_physical(n) - Z_physical(n+1))];
end
a_physical(i:i:end)=[]; 
a_virtual(i:i:end)=[]; 

a_physical_vector(i:i:end,:)=[];
a_virtual_vector(i:i:end,:)=[];

for n=1:(i.*(j-1))
    d_physical(n)=sqrt((X_physical(n) - X_physical(n+i)).^2 +(Y_physical(n)-Y_physical(n+i)).^2+(Z_physical(n)-Z_physical(n+i)).^2);
    d_virtual(n)=sqrt((X_virtual(n) - X_virtual(n+i)).^2 +(Y_virtual(n)-Y_virtual(n+i)).^2+(Z_virtual(n)-Z_virtual(n+i)).^2);
    d_virtual_vector(n,:)=[(X_virtual(n)-X_virtual(n+i)),(Y_virtual(n)-Y_virtual(n+i)),(Z_virtual(n)-Z_virtual(n+i))];
    d_physical_vector(n,:)=[(X_physical(n)-X_physical(n+i)),(Y_physical(n)-Y_physical(n+i)),(Z_physical(n)-Z_physical(n+i))];
end
i_physical = a_physical;
i_physical_vector=a_physical_vector;
j_physical =d_physical;
j_physical_vector=d_physical_vector;
i_virtual = a_virtual;
i_virtual_vector=a_virtual_vector;
j_virtual =d_virtual;
j_virtual_vector=d_virtual_vector;

c_physical = a_physical;
c_physical_vector=a_physical_vector;
c_virtual = a_virtual; 
c_virtual_vector=a_virtual_vector;
b_physical = d_physical; 
b_physical_vector=d_physical_vector;
b_virtual = d_virtual;
b_virtual_vector=d_virtual_vector;

a_physical(((i-1).*(j-1))+1:((i-1).*j))=[];
a_physical_vector(((i-1).*(j-1))+1:((i-1).*j),:)=[];
a_virtual(((i-1).*(j-1))+1:((i-1).*j))=[];
a_virtual_vector(((i-1).*(j-1))+1:((i-1).*j),:)=[];
c_physical(1:(i-1))=[];
c_physical_vector(1:(i-1),:)=[];
c_virtual(1:(i-1))=[];
c_virtual_vector(1:(i-1),:)=[];
b_physical(1:i:end)=[];
b_physical_vector(1:i:end,:)=[];
b_virtual(1:i:end)=[];
b_virtual_vector(1:i:end,:)=[];
d_physical(i:i:end)=[]; 
d_physical_vector(i:i:end,:)=[];
d_virtual(i:i:end)=[];
d_virtual_vector(i:i:end,:)=[];


% %note that the cross product is going to depend on intial mesh orientation
n_physical=cross(a_physical_vector,d_physical_vector);
n_virtual=cross(a_virtual_vector,d_virtual_vector);


for n=1:length(n_physical)
%     n11_physical(n)=n_physical(n,1).*n_physical(n,1);
%     n12_physical(n)=n_physical(n,1).*n_physical(n,2);
%     n13_physical(n)=n_physical(n,1).*n_physical(n,3);
%     n21_physical(n)=n_physical(n,2).*n_physical(n,1);
%     n22_physical(n)=n_physical(n,2).*n_physical(n,2);
%     n23_physical(n)=n_physical(n,2).*n_physical(n,3);
%     n31_physical(n)=n_physical(n,3).*n_physical(n,1);
%     n32_physical(n)=n_physical(n,3).*n_physical(n,2);
%     n33_physical(n)=n_physical(n,3).*n_physical(n,3);
%     n11_virtual(n)=n_virtual(n,1).*n_virtual(n,1);
%     n12_virtual(n)=n_virtual(n,1).*n_virtual(n,2);
%     n13_virtual(n)=n_virtual(n,1).*n_virtual(n,3);
%     n21_virtual(n)=n_virtual(n,2).*n_virtual(n,1);
%     n22_virtual(n)=n_virtual(n,2).*n_virtual(n,2);
%     n23_virtual(n)=n_virtual(n,2).*n_virtual(n,3);
%     n31_virtual(n)=n_virtual(n,3).*n_virtual(n,1);
%     n32_virtual(n)=n_virtual(n,3).*n_virtual(n,2);
%     n33_virtual(n)=n_virtual(n,3).*n_virtual(n,3);
    norm_physical(n)=norm(n_physical(n,:));
    norm_virtual(n)=norm(n_virtual(n,:));
    cap_n_physical(n,:)=n_physical(n,:)./(norm_physical(n));
    cap_n_virtual(n,:)=n_virtual(n,:)./(norm_virtual(n));
end

%cap_n_physical=n_physical./(norm_physical);
%cap_n_virtual=n_virtual./(norm_virtual);

% norm_physical_cell=cell(i-1,j-1);
% norm_virtual_cell=cell(i-1,j-1);
% 
% %This Cell Population method is incorrect
%   for n=1:(j-1)
%     for m=1:(i-1)
%         norm_physical_cell{n,m}=norm_physical(n.*m).*[n11_physical(n.*m),n12_physical(n.*m),n13_physical(n.*m);n21_physical(n.*m),n22_physical(n.*m),n23_physical(n.*m);n31_physical(n.*m),n32_physical(n.*m),n33_physical(n.*m)];
%         norm_virtual_cell{n,m}=norm_virtual(n.*m).*[n11_virtual(n.*m),n12_virtual(n.*m),n13_virtual(n.*m);n21_virtual(n.*m),n22_virtual(n.*m),n23_virtual(n.*m);n31_virtual(n.*m),n32_virtual(n.*m),n33_virtual(n.*m)];
%     end
%   end
% 
% 
% euclidean_metric_cell=cell(i-1,j-1);
% for n=1:j-1
%     for m=1:i-1
%         euclidean_metric_cell{n,m}=[1,0,0;0,1,0;0,0,1];
%     end
% end
% 
%  for n=1:(j-1)
%     for m=1:(i-1)
%         h_physical_cell{n,m}=euclidean_metric_cell{n,m}-norm_physical_cell{n,m};
%         h_virtual_cell{n,m}=euclidean_metric_cell{n,m}-norm_virtual_cell{n,m};
%     end
%   end
% 
% 
% edge_length_physical = [a_physical;b_physical;c_physical;d_physical]';
% edge_length_virtual = [a_virtual;b_virtual;c_virtual;d_virtual]';
% P_physical = a_physical+b_physical+c_physical+d_physical;
% P_virtual = a_virtual+b_virtual+c_virtual+d_physical;

for n = 1:(i.*j)-1
    h11(n) = 1 + ((Z_virtual(n)-Z_virtual(n+1))./(X_virtual(n)-X_virtual(n+1))).^2;
    h12_x(n) = (Z_virtual(n)-Z_virtual(n+1))./(X_virtual(n)-X_virtual(n+1));
end
h11(i:i:end)=[]; h11(((i-1).*(j-1))+1:((i-1).*j))=[];
h12_x(i:i:end)=[]; h12_x(((i-1).*(j-1))+1:((i-1).*j))=[];
for n = 1:(i.*(j-1))
    h22(n) = 1 + ((Z_virtual(n)-Z_virtual(n+i))./(Y_virtual(n)-Y_virtual(n+i))).^2;
    h12_y(n)= (Z_virtual(n)-Z_virtual(n+i))./(Y_virtual(n)-Y_virtual(n+i));
end
h22(i:i:end)=[];
h12_y(i:i:end)=[];

h12 = h12_x.*h12_y;

%new population technique that assumes that x and y are not independent
for n = 1:(i.*j)-1
    h11_new(n) = 1 + ((Y_virtual(n)-Y_virtual(n+1))./(X_virtual(n)-X_virtual(n+1))).^2 + ((Z_virtual(n)-Z_virtual(n+1))./(X_virtual(n)-X_virtual(n+1))).^2;
    h12_C_x_new(n) = (Z_virtual(n)-Z_virtual(n+1))./(X_virtual(n)-X_virtual(n+1));
    h12_B_x_new(n) = (Y_virtual(n)-Y_virtual(n+1))./(X_virtual(n)-X_virtual(n+1));
    lambda11(n)=(X_virtual(n)-X_virtual(n+1))./(X_physical(n)-X_physical(n+1));
    lambda12(n)=(X_virtual(n)-X_virtual(n+1))./(Y_physical(n)-Y_physical(n+1));
end
h11_new(i:i:end)=[]; h11_new(((i-1).*(j-1))+1:((i-1).*j))=[];
lambda11(i:i:end)=[];lambda11(((i-1).*(j-1))+1:((i-1).*j))=[];
lambda12(i:i:end)=[];lambda12(((i-1).*(j-1))+1:((i-1).*j))=[];
h12_C_x_new(i:i:end)=[]; h12_C_x_new(((i-1).*(j-1))+1:((i-1).*j))=[];
h12_B_x_new(i:i:end)=[]; h12_B_x_new(((i-1).*(j-1))+1:((i-1).*j))=[];
for n = 1:(i.*(j-1))
    h22_new(n) = ((X_virtual(n)-X_virtual(n+i))./(Y_virtual(n)-Y_virtual(n+i))).^2 + 1 + ((Z_virtual(n)-Z_virtual(n+i))./(Y_virtual(n)-Y_virtual(n+i))).^2;
    h12_C_y_new(n)= (Z_virtual(n)-Z_virtual(n+i))./(Y_virtual(n)-Y_virtual(n+i));
    h12_A_y_new(n) = (X_virtual(n)-X_virtual(n+i))./(Y_virtual(n)-Y_virtual(n+i));
    lambda21(n)=(Y_virtual(n)-Y_virtual(n+i))./(X_physical(n)-X_physical(n+i));
    lambda22(n)=(Y_virtual(n)-Y_virtual(n+i))./(Y_physical(n)-Y_physical(n+i));
end
h22_new(i:i:end)=[];
h12_C_y_new(i:i:end)=[];
h12_A_y_new(i:i:end)=[];
lambda21(i:i:end)=[];
lambda22(i:i:end)=[];

h12_C = h12_C_x_new.*h12_C_y_new;
h12_A = h12_A_y_new; h12_B = h12_B_x_new;
h12_new = h12_A + h12_B + h12_C;

for n=1:cells
    E_prime(n) = (1./((h11_new(n).*h22_new(n) - (h12(n)).^2))).*h22_new(n);
    F_prime(n)= (1./((h11_new(n).*h22_new(n) - (h12(n)).^2))).*-.1*h12(n);
    G_prime(n)= (1./((h11_new(n).*h22_new(n) - (h12(n)).^2))).*h11_new(n);
end


%Check to see about i and j
for n=1:(i.*j-i)-1
    diag2_physical_length(n)=sqrt((X_physical(n)-X_physical(n+i+1)).^2 + (Y_physical(n)-Y_physical(n+i+1)).^2 + (Z_physical(n)-Z_physical(n+i+1)).^2);
    diag2_virtual_length(n)=sqrt((X_virtual(n)-X_virtual(n+i+1)).^2 + (Y_virtual(n)-Y_virtual(n+i+1)).^2 + (Z_virtual(n)-Z_virtual(n+i+1)).^2);
end
diag2_physical_length(i:i:end)=[];diag2_virtual_length(i:i:end)=[];

for n=1:(i.*j-i)
    diag1_physical_length(n)=sqrt((X_physical(n+1)-X_physical(n+i)).^2 + (Y_physical(n+1)-Y_physical(n+i)).^2 + (Z_physical(n+1)-Z_physical(n+i)).^2);
    diag1_virtual_length(n)=sqrt((X_virtual(n+1)-X_virtual(n+i)).^2 +(Y_virtual(n+1)-Y_virtual(n+i)).^2 + (Z_virtual(n+1)-Z_virtual(n+i)).^2);
end
diag1_physical_length(i:i:end)=[];diag1_virtual_length(i:i:end)=[];

diag_length_physical = [diag1_physical_length;diag2_physical_length]';
diag_length_virtual = [diag1_virtual_length;diag2_virtual_length]';
%Interior angle calculations

%alpha:corner of a and d
%beta:corner of a and b
%gamma:corner of b and c
%delta:corner of c and d
 for n=1:(i-1)*(j-1)
    alpha_physical(n)=2.*atand(sqrt(((diag1_physical_length(n)).^2-(d_physical(n)-a_physical(n)).^2)./((d_physical(n)+a_physical(n)).^2 - (diag1_physical_length(n)).^2)));
    alpha_virtual(n)=2.*atand(sqrt(((diag1_virtual_length(n)).^2-(d_virtual(n)-a_virtual(n)).^2)./((d_virtual(n)+a_virtual(n)).^2 - (diag1_virtual_length(n)).^2)));
    beta_physical(n)=2.*atand(sqrt(((diag2_physical_length(n)).^2-(a_physical(n)-b_physical(n)).^2)./((a_physical(n)+b_physical(n)).^2 - (diag2_physical_length(n)).^2)));
    beta_virtual(n)=2.*atand(sqrt(((diag2_virtual_length(n)).^2-(a_virtual(n)-b_virtual(n)).^2)./((a_virtual(n)+b_virtual(n)).^2 - (diag2_virtual_length(n)).^2)));
    gamma_physical(n)=2.*atand(sqrt(((diag1_physical_length(n)).^2-(b_physical(n)-c_physical(n)).^2)./((b_physical(n)+c_physical(n)).^2 - (diag1_physical_length(n)).^2)));
    gamma_virtual(n)=2.*atand(sqrt(((diag1_virtual_length(n)).^2-(b_virtual(n)-c_virtual(n)).^2)./((b_virtual(n)+c_virtual(n)).^2 - (diag1_virtual_length(n)).^2)));
    delta_physical(n)=2.*atand(sqrt(((diag2_physical_length(n)).^2-(c_physical(n)-d_physical(n)).^2)./((c_physical(n)+d_physical(n)).^2 - (diag2_physical_length(n)).^2)));
    delta_virtual(n)=2.*atand(sqrt(((diag2_virtual_length(n)).^2-(c_virtual(n)-d_virtual(n)).^2)./((c_virtual(n)+d_virtual(n)).^2 - (diag2_virtual_length(n)).^2)));
 end

 %%Geometric Method
 %This particular definition of the the metric tensor relies on edge
 %lengths a and d, as well as the interior angle alpha
 
for n=1:length(a_physical_vector)
    E_geometric(n) = dot(a_physical_vector(n),a_physical_vector(n)); G_geometric(n) = dot(d_physical_vector(n),d_physical_vector(n));
    F_geometric(n) = dot(a_physical_vector(n),d_physical_vector(n));
end

for n=1:length(a_physical)
    E_geometric_other(n)=(a_physical(n)).^2; G_geometric_other(n)=(d_physical(n)).^2;
    F_geometric_other(n)=(a_physical(n)).*(d_physical(n)).*cosd(alpha_physical(n));
end

angles_physical = [alpha_physical;beta_physical;gamma_physical;delta_physical]';
angles_virtual = [alpha_virtual;beta_virtual;gamma_virtual;delta_virtual]';

physical_sum=(alpha_physical+beta_physical+gamma_physical+delta_physical)';
virtual_sum=(alpha_virtual+beta_virtual+gamma_virtual+delta_virtual)';

%Tri1=Triangle(bcdiagl) Tri2=Triangle(addiag1)
% k=(a+b+c)./2;
for n=1:length(diag1_physical_length)
    k_Tri1_physical(n)=(b_physical(n)+c_physical(n)+diag1_physical_length(n))./2;
    k_Tri1_virtual(n)=(b_virtual(n)+c_virtual(n)+diag1_virtual_length(n))./2;
    k_Tri2_physical(n)=(a_physical(n)+d_physical(n)+diag1_physical_length(n))./2;
    k_Tri2_virtual(n)=(a_virtual(n)+d_virtual(n)+diag1_virtual_length(n))./2;
end

%Heron's Formula: A=sqrt((k(k-a)(k-b)(k-c))

for n=1:length(diag1_physical_length)
    A1_physical_diag1(n)=sqrt((k_Tri1_physical(n).*(k_Tri1_physical(n)-b_physical(n)).*(k_Tri1_physical(n)-c_physical(n)).*(k_Tri1_physical(n)-diag1_physical_length(n))));
    A1_virtual_diag1(n)=sqrt((k_Tri1_virtual(n).*(k_Tri1_virtual(n)-b_virtual(n)).*(k_Tri1_virtual(n)-c_virtual(n)).*(k_Tri1_virtual(n)-diag1_virtual_length(n))));
    A2_physical_diag1(n)=sqrt((k_Tri2_physical(n).*(k_Tri2_physical(n)-a_physical(n)).*(k_Tri2_physical(n)-d_physical(n)).*(k_Tri2_physical(n)-diag1_physical_length(n))));
    A2_virtual_diag1(n)=sqrt((k_Tri2_virtual(n).*(k_Tri2_virtual(n)-a_virtual(n)).*(k_Tri2_virtual(n)-d_virtual(n)).*(k_Tri2_virtual(n)-diag1_virtual_length(n))));
end

A_physical_diag1=(A1_physical_diag1+A2_physical_diag1); A_virtual_diag1=A1_virtual_diag1+A2_virtual_diag1;


%Test againts diag2
%Tri1=Triangle(abdiag2) Tri2=Triangle(cddiag2)
for n=1:length(diag1_physical_length)
    k_Tri1_physical_diag2(n)=(a_physical(n)+b_physical(n)+diag2_physical_length(n))./2;
    k_Tri1_virtual_diag2(n)=(a_virtual(n)+b_virtual(n)+diag2_virtual_length(n))./2;
    k_Tri2_physical_diag2(n)=(c_physical(n)+d_physical(n)+diag2_physical_length(n))./2;
    k_Tri2_virtual_diag2(n)=(c_virtual(n)+d_virtual(n)+diag2_virtual_length(n))./2;
end

%Heron's Formula: A=sqrt((k(k-a)(k-b)(k-c))

for n=1:length(diag2_physical_length)
    A1_physical_diag2(n)=sqrt((k_Tri1_physical_diag2(n).*(k_Tri1_physical_diag2(n)-a_physical(n)).*(k_Tri1_physical_diag2(n)-b_physical(n)).*(k_Tri1_physical_diag2(n)-diag2_physical_length(n))));
    A1_virtual_diag2(n)=sqrt((k_Tri1_virtual_diag2(n).*(k_Tri1_virtual_diag2(n)-a_virtual(n)).*(k_Tri1_virtual_diag2(n)-b_virtual(n)).*(k_Tri1_virtual_diag2(n)-diag2_virtual_length(n))));
    A2_physical_diag2(n)=sqrt((k_Tri2_physical_diag2(n).*(k_Tri2_physical_diag2(n)-c_physical(n)).*(k_Tri2_physical_diag2(n)-d_physical(n)).*(k_Tri2_physical_diag2(n)-diag2_physical_length(n))));
    A2_virtual_diag2(n)=sqrt((k_Tri2_virtual_diag2(n).*(k_Tri2_virtual_diag2(n)-c_virtual(n)).*(k_Tri2_virtual_diag2(n)-d_virtual(n)).*(k_Tri2_virtual_diag2(n)-diag2_virtual_length(n))));
end

%Calculation of Arbitrary triangle
Alpha_diag2_1=sqrt((a_physical+b_physical).^2 - diag2_physical_length);
Beta_diag2_1=sqrt((diag2_physical_length).^2 - (a_physical-b_physical).^2);
Alpha_diag2_2=sqrt((c_physical+d_physical).^2 - diag2_physical_length);
Beta_diag2_2=sqrt((diag2_physical_length).^2 - (c_physical-d_physical).^2);

for n=1:length(diag2_physical_length)
    A1_physical_diag2_other(n)=(Alpha_diag2_1(n).*Beta_diag2_1(n))./4;
    A2_physical_diag2_other(n)=(Alpha_diag2_2(n).*Beta_diag2_2(n))./4;
end

A_physical_diag2_other = A1_physical_diag2_other+A2_physical_diag2_other;

A_physical_diag2=A1_physical_diag2+A2_physical_diag2; A_virtual_diag2=A1_virtual_diag2+A2_virtual_diag2;

%Calculating Area using General Quadrilateral Approach
p_physical = a_physical + b_physical + c_physical + d_physical;
p_virtual = a_virtual + b_virtual + c_virtual + d_virtual;
k_physical = p_physical./2; k_virtual = p_virtual./2;
phi_physical = (alpha_virtual+gamma_virtual)./2; phi_virtual = (alpha_virtual+gamma_virtual)./2;
%%A_gen_quad = sqrt((k-a).*(k-b).*(k-c).*(k-d)-(a.*b.*c.*d).*((cosd(phi)).^2));
abcd_physical=(a_physical.*b_physical.*c_physical.*d_physical); abcd_virtual = (a_virtual.*b_virtual.*c_virtual.*d_virtual);
ka_physical = k_physical - a_physical; kb_physical = k_physical-b_physical; kc_physical=k_physical-c_physical; kd_physical=k_physical-d_physical;
ka_virtual = k_virtual - a_virtual; kb_virtual = k_virtual - b_virtual; kc_virtual = k_virtual - c_virtual; kd_virtual = k_virtual- d_virtual;

A_quad_physical = sqrt(((ka_physical).*(kb_physical).*(kc_physical).*(kd_physical))-(abcd_physical.*(cosd(phi_physical)).^2));
A_quad_virtual = sqrt(((ka_virtual).*(kb_virtual).*(kc_virtual).*(kd_virtual))-(abcd_virtual.*(cosd(phi_virtual)).^2));

% A_physical_perc_diff=(abs((A_physical_diag1-A_physical_diag2)./((A_physical_diag1+A_physical_diag2)./2))).*100;
% A_virtual_perc_diff=(abs((A_virtual_diag1-A_virtual_diag2)./((A_virtual_diag1+A_virtual_diag2)./2))).*100;

%Reshaping A_quad_physical and A_quad_virtual
% figure
% A_quad_ratio = eps_virtual.*(A_quad_virtual./A_quad_physical);
% A_quad_ratio_matrix = reshape(A_quad_ratio,[],j-1)';
% imagesc(flipud(A_quad_ratio_matrix));
% title(sprintf('Permittivity Map (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along i-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,i-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along j-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,j-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar
% 
% figure
% A_physical_diag1_matrix = reshape(A_physical_diag1,[],j-1)';
% imagesc(flipud(A_physical_diag1_matrix));
% title(sprintf('Physical Space Area Distribution',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along i-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,i-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along j-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,j-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar
% 
% figure
% A_virtual_diag1_matrix = reshape(A_virtual_diag1,[],j-1)';
% imagesc(flipud(A_virtual_diag1_matrix));
% title(sprintf('Virtual Space Area Distribution',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along i-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,i-1,11),'XTickLabel',(linspace(min(X_virtual),max(X_virtual),11)),'Fontweight','bold')
% ylabel('Pos. along j-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,j-1,11),'YTickLabel',(linspace(max(Y_virtual),min(Y_virtual),11)))
% colorbar

% for n=1:(i-1).*(j-1)
%     physical_perm_matrix_true(n,:)=physical_perm_matrix_other(n,:).*Virtual_perm(n,:);
% end
% 
% physical_perm_eps_11=physical_perm_matrix_true(:,1);
% 
% fid=fopen('true_epsilon.txt','w+');
% for n=1:(i-1).*(j-1)
%     for m=1:4
%         fprintf(fid,'%f',physical_perm_matrix_true(n,m));
%         fprintf(fid,'\t');
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);

%Metric tensor for each cell
%g1=a; g2=d; theta = alpha
% g1_physical = a_physical;
% g1_physical_vector=a_physical_vector;
% g1_virtual=a_virtual;
% g1_virtual_vector=a_virtual_vector;
% g2_physical = d_physical;
% g2_physical_vector=d_physical_vector;
% g2_virtual=d_virtual;
% g2_virtual_vector=d_virtual_vector;
% 
% 
% for n=1:length(g1_physical_vector)
%     g11_physical(n)=dot(g1_physical_vector(n,:),g1_physical_vector(n,:));
%     g22_physical(n)=dot(g2_physical_vector(n,:),g2_physical_vector(n,:));
%     g12_physical(n)=dot(g1_physical_vector(n,:),g2_physical_vector(n,:));
%     g21_physical(n)=dot(g2_physical_vector(n,:),g1_physical_vector(n,:));
%     g11_virtual(n)=dot(g1_virtual_vector(n,:),g1_virtual_vector(n,:));    
%     g22_virtual(n)=dot(g2_virtual_vector(n,:),g2_virtual_vector(n,:));
%     g12_virtual(n)=dot(g1_virtual_vector(n,:),g2_virtual_vector(n,:));
%     g21_virtual(n)=dot(g2_virtual_vector(n,:),g1_virtual_vector(n,:));
% 
% end

% physical_metric_area=sqrt((g11_physical_vector).*(g22_physical_vector)+((g12_physical_vector).*(g21_physical_vector)));
% virtual_metric_area=sqrt((g11_virtual_vector).*(g22_virtual_vector)+((g12_virtual_vector).*(g21_virtual_vector)));

%Keep in mind that this equation is only valid for parallelograms and not
%arbitrart quadrilaterals

% physical_metric_area_revised=sqrt((g11_physical.*g22_physical)-(g12_physical).^2);
% virtual_metric_area_revised=sqrt((g11_virtual.*g22_virtual)-(g12_virtual).^2);
% 
% metric_alpha_physical=acosd(g12_physical./(sqrt(g11_physical.*g22_physical)));
% metric_alpha_virtual=acosd(g12_virtual./(sqrt(g11_virtual.*g22_virtual)));
% 
% % physical_percent_difference = abs(physical_metric_area-A_quad_physical)./((abs(physical_metric_area+A_quad_physical)./2));
% % virtual_percent_difference = abs(virtual_metric_area-A_quad_virtual)./((abs(virtual_metric_area+A_quad_virtual)./2));
% 
% physical_percent_difference_new=abs(physical_metric_area_revised-A_quad_physical)./(abs(physical_metric_area_revised+A_quad_physical)./2);
% virtual_percent_difference_new=abs(virtual_metric_area_revised-A_quad_virtual)./(abs(virtual_metric_area_revised+A_quad_virtual)./2);
% 
% 
% g_physical=cell((i-1),(j-1));
% g_virtual=cell((i-1),(j-1));
% g_physical_vector=cell((i-1),(j-1));
% g_virtual_vector=cell((i-1),(j-1));

%Redo this population method
% for n=1:(j-1)
%     for m=1:(i-1)
%         g_physical{n,m}=[g11_physical(n.*m),g12_physical(n.*m);g21_physical(n.*m),g22_physical(n.*m)];
%         g_virtual{n,m}=[g11_virtual(n.*m),g12_virtual(n.*m);g21_virtual(n.*m),g22_virtual(n.*m)];
%         g_physical_vector{n,m}=[g11_physical_vector(n.*m),g12_physical_vector(n.*m);g21_physical_vector(n.*m),g22_physical_vector(n.*m)];
%         g_virtual_vector{n,m}=[g11_virtual_vector(n.*m),g12_virtual_vector(n.*m);g21_virtual_vector(n.*m),g22_virtual_vector(n.*m)];
%     end
% end
% 
% for n=1:(j-1)
%     for m=1:(i-1)
%         g_physical_determinant(n,m)=det(g_physical{n,m});
%         g_virtual_determinant(n,m)=det(g_virtual{n,m});
%         g_physical_determinant_vector(n,m)=det(g_physical_vector{n,m});
%         g_virtual_determinant_vector(n,m)=det(g_virtual_vector{n,m});
%     end
% end
% 
% physical_perm = cell((i-1),(j-1));
% for n=1:(j-1)
%     for m=1:(i-1)
%         physical_perm{n,m}=((sqrt(g_physical_determinant(n,m)))./(sqrt(g_virtual_determinant(n,m)))).*g_physical{n,m};
%         physical_perm_other{n,m}=(A_quad_physical(n.*m)./A_quad_virtual(n.*m).*g_physical{n,m});
%     end
% end

%.*g_physical{n,m}

% for n=1:(i-1).*(j-1)
%     physical_perm_matrix(n,:)=reshape(cell2mat(physical_perm(n)),1,[]);
%     physical_perm_matrix_other(n,:)=reshape(cell2mat(physical_perm_other(n)),1,[]); 
% end


% physical_perm_eps_11=physical_perm_matrix_other(:,1);
% physical_perm_eps_12=physical_perm_matrix_other(:,2);
% physical_perm_eps_21=physical_perm_matrix_other(:,3);
% physical_perm_eps_22=physical_perm_matrix_other(:,4);
% 
% info_physical=cell((i-1),(j-1));
% info_virtual=cell((i-1),(j-1));

% for n=1:(j-1)
%     for m=1:(i-1)
%         
%         info_physical{n,m}=[a_physical(n.*m),b_physical(n.*m),c_physical(n.*m),d_physical(n.*m); alpha_physical(n.*m), beta_physical(n.*m), gamma_physical(n.*m),delta_physical(n.*m);diag1_physical_length(n.*m),diag2_physical_length(n.*m),A_physical_diag1(n.*m),physical_sum(n.*m)];
%         info_virtual{n,m}=[a_virtual(n.*m),b_virtual(n.*m),c_virtual(n.*m),d_virtual(n.*m); alpha_virtual(n.*m), beta_virtual(n.*m), gamma_virtual(n.*m),delta_virtual(n.*m);diag1_virtual_length(n.*m),diag2_virtual_length(n.*m),A_virtual_diag1(n.*m),virtual_sum(n.*m)];
% 
%     end
% end


%g_physical=[g11_physical,g12_physical,g21_physical,g22_physical];
%g_virtual =[g11_virtual,g12_virtual; g21_virtual,g22_virtual];
%C=mat2cell(g_physical,[2 51200], [2 2])

% x=1:length(A_physical_perc_diff);
% figure
% subplot(3,2,1)
% scatter3(X_virtual,Y_virtual,Z_virtual,15,'rd')
% title(sprintf('Virtual Space (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Pos. along X-axis (mm)','FontSize',11,'Fontweight','bold')
% ylabel('Pos. along Y-axis (mm)','FontSize',11,'Fontweight','bold')
% zlabel('Pos. along Z-axis (mm)','FontSize',11,'Fontweight','bold')
% 
% subplot(3,2,2)
% scatter3(X_physical,Y_physical,Z_physical,15,'bd')
% title(sprintf('Physical Space (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Pos. along X-axis (mm)','FontSize',11,'Fontweight','bold')
% ylabel('Pos. along Y-axis (mm)','FontSize',11,'Fontweight','bold')
% zlabel('Pos. along Z-axis (mm)','FontSize',11,'Fontweight','bold')
% 
% subplot(3,2,3)
% plot(x,A_virtual_perc_diff,'r')
% title(sprintf('Virtual Space (%dx%d Surface Mesh)',i-1,j-1),'FontSize',12,'Fontweight','bold')
% xlabel('Element index','FontSize',11,'Fontweight','bold')
% ylabel('Percentage Difference (A_D_1,A_D_2)','FontSize',11,'Fontweight','bold')
% xlim([1 length(virtual_sum)])
% 
% subplot(3,2,4)
% plot(x,A_physical_perc_diff,'b')
% title(sprintf('Physical Space (%dx%d Surface Mesh)',i-1,j-1),'FontSize',12,'Fontweight','bold')
% xlabel('Element index','FontSize',11,'Fontweight','bold')
% ylabel('Percentage Difference (A_D_1,A_D_2)','FontSize',11,'Fontweight','bold')
% xlim([1 length(physical_sum)])
% 
% A_physical_perc_diff_matrix=reshape(A_physical_perc_diff,[],j-1);
% A_virtual_perc_diff_matrix=reshape(A_virtual_perc_diff,[],j-1);


% 
% subplot(3,2,5)
% plot(x,virtual_sum,'r')
% title(sprintf('Summation of interior angles in Virtual Space for %dx%d Surface Mesh',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Element index','FontSize',11,'Fontweight','bold')
% ylabel('Degrees','FontSize',11,'Fontweight','bold')
% xlim([1 length(virtual_sum)])
% 
% subplot(3,2,6)
% plot(x,physical_sum,'b')
% title(sprintf('Summation of interior angles in Physical Space for %dx%d Surface Mesh',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Element index','FontSize',11,'Fontweight','bold')
% ylabel('Degrees','FontSize',11,'Fontweight','bold')
% xlim([1 length(physical_sum)])

A = h11_new; B=h12_new; C=h22_new; W=lambda11; X=lambda12; Y=lambda21; Z=lambda22;
denom = W.*Z - X.*Y;

% for n = length(h11_new)
%     gamma11(n)=X(n).*((A(n).*X(n))-
%     
physical_perm_eps_11 = ((A_virtual_diag1)./(A_physical_diag1)).*(E_prime);
physical_perm_eps_12 = ((A_virtual_diag1)./(A_physical_diag1)).*(F_prime);
physical_perm_eps_22 = ((A_virtual_diag1)./(A_physical_diag1)).*(G_prime);

 eps_11_matrix=reshape(physical_perm_eps_11,[],j-1);
 eps_12_matrix=reshape(physical_perm_eps_12,[],j-1);
% eps_21_matrix=reshape(physical_perm_eps_21,[],j-1);
 eps_22_matrix=reshape(physical_perm_eps_22,[],j-1);
 
%  perm_master=A_quad_physical./A_quad_virtual;
%  perm_master=reshape(perm_master,[],j-1);

% figure
% imagesc(flipud(perm_master)')
% title(sprintf('\\epsilon_1 master (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along i-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along j-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar

%caxis([1 5])
figure
imagesc(flipud(eps_11_matrix'))
title(sprintf('\\epsilon_1_1, \\mu_1_1 (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
xlim([1 i-1])
ylim([1 j-1])
xlabel('X','FontSize',11,'Fontweight','bold')
set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
ylabel('Y','FontSize',11,'Fontweight','bold')
set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
colorbar
% 
figure
imagesc(flipud(eps_12_matrix'))
title(sprintf('\\epsilon_1_2, \\mu_1_2 (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
xlim([1 i-1])
ylim([1 j-1])
xlabel('X','FontSize',11,'Fontweight','bold')
set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
ylabel('Y','FontSize',11,'Fontweight','bold')
set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
colorbar

% figure
% imagesc(flipud(eps_21_matrix))
% title(sprintf('\\epsilon_2_1 (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along i-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along j-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar

E_prime_new_matrix =reshape(h11_new,[],j-1);
F_prime_new_matrix=reshape(h12_new,[],j-1);
G_prime_new_matrix=reshape(h22_new,[],j-1);

% figure
% imagesc(flipud(E_prime_new_matrix'))
% title(sprintf('E (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along x-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along y-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar
% grid on
% set(gcf, 'Color', 'w')

% figure
% imagesc(flipud(F_prime_new_matrix'))
% title(sprintf('F (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along x-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along y-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar
% grid on
% set(gcf, 'Color', 'w')

% figure
% imagesc(flipud(G_prime_new_matrix'))
% title(sprintf('G (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along x-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along y-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar
% grid on
% set(gcf, 'Color', 'w')

figure
imagesc(flipud(eps_22_matrix'))
title(sprintf('\\epsilon_2_2, \\mu_2_2 (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
xlim([1 i-1])
ylim([1 j-1])
xlabel('X','FontSize',11,'Fontweight','bold')
set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
ylabel('Y','FontSize',11,'Fontweight','bold')
set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
colorbar

% x=1:length(physical_sum);
% figure
% plot(x,physical_sum,'b','Linewidth',2)
% title(sprintf('Summation of Interior Angles (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Patch index','FontSize',11,'Fontweight','bold')
% ylabel('degrees','FontSize',11,'Fontweight','bold')
% xlim([1 length(physical_sum)])
% 
% summation_of_angles = figure();
% set(summation_of_angles,'Position',[100 100 600 400]);
% plot(x,physical_sum,'k','LineWidth',3)
% title(sprintf('Summation of interior angles (%dx%d Surface Mesh)',i-1,j-1),'FontSize',18,'Interpreter','Latex')
% xlabel('Patch index','FontSize',16,'Interpreter','Latex','FontWeight','bold')
% ylabel('Degrees','FontSize',16,'Fontweight','bold','Interpreter','Latex')
% xlabel([1 length(physical_sum)])
% grid on
% set(gcf, 'Color', 'w')



% figure
% plot(x,virtual_sum,'r')
% title(sprintf('Virtual Space for %dx%d Surface Mesh',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Element index','FontSize',11,'Fontweight','bold')
% ylabel('Summation of Interior Angles (Degrees)','FontSize',11,'Fontweight','bold')
% xlim([1 length(virtual_sum)])
% 
% figure
% plot(x,physical_sum,'b')
% title(sprintf('Physical Space for %dx%d Surface Mesh',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Element index','FontSize',11,'Fontweight','bold')
% ylabel('Summation of Interior Angles (Degrees)','FontSize',11,'Fontweight','bold')
% xlim([1 length(physical_sum)])


% ref_index_map = sqrt(A_virtual_diag1./A_physical_diag1);
% ref_map_matrix=reshape(ref_index_map,[],j-1)';
% 
% mag_virtual_matrix = reshape(A_virtual_diag1,i-1,j-1);
% mag_physical_matrix = reshape(A_physical_diag1,i-1,j-1)';
% 
% Area_ratio= mag_virtual_matrix./mag_physical_matrix;
% Area_ratio_line = A_virtual_diag1./A_physical_diag1;

% figure
% hist(Area_ratio_line,100)
% title(sprintf('Virtual Area/Physical Area Frequency(%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Area ratio','FontSize',11,'Fontweight','bold')
% ylabel('Frequency','FontSize',11,'Fontweight','bold')

% figure
% plot(x,Area_ratio_line,'b')
% title(sprintf('Virtual Area/Physical Area (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Cell index','FontSize',11,'Fontweight','bold')
% ylabel('Area ratio','FontSize',11,'Fontweight','bold')
% xlim([1 length(physical_sum)])
% 
% figure
% imagesc(flipud(Area_ratio))
% title(sprintf('Virtual Area/Physical Area (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along i-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along j-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar
% 
% figure
% imagesc(flipud(mag_virtual_matrix))
% title(sprintf('Virtual Space Area Map for %dx%d Surface Mesh',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along i-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along j-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar
% 
% figure
% imagesc(flipud(mag_physical_matrix))
% title(sprintf('Physical Space Area Map for %dx%d Surface Mesh',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along i-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along j-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar 


% subplot(3,3,8)
% imagesc(flipud(ref_map_matrix))
% title('Refractive Index of Physical Space (Diag1)','FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along X-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along Y-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar 
% 
% ref_map_diag2 = sqrt(A_virtual_diag2./A_physical_diag2);
% ref_map_matrix_diag2=reshape(ref_map_diag2,[],j-1)';
% 
% subplot(3,3,9)
% imagesc(flipud(ref_map_matrix_diag2))
% title('Refractive Index of Physical Space (Diag2)','FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along X-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_physical),max(X_physical),11)),'Fontweight','bold')
% ylabel('Pos. along Y-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_physical),min(Y_physical),11)))
% colorbar 

% fid = fopen('epsilon.txt','w');
% for n=1:j-1
%     for m=1:i-1
%         
%     fprintf(fid, '%f', mag_virtual_matrix(n,m));
%     fprintf(fid,'\t');
%     end
%     
%     fprintf(fid, '\n');
% end


% PAY ATTENTION TO THIS
eps_11_matrix = eps_11_matrix';
eps_12_matrix = eps_12_matrix';
eps_22_matrix = eps_22_matrix';
fid = fopen('eps_mu_11.txt','w');
for n=1:j-1
    for m=1:i-1
        
    fprintf(fid, '%f', eps_11_matrix(n,m));
    fprintf(fid,'\t');
    end
    
    fprintf(fid, '\n');
end

fid = fopen('eps_mu_12.txt','w');
for n=1:j-1
    for m=1:i-1
        
    fprintf(fid, '%f', eps_12_matrix(n,m));
    fprintf(fid,'\t');
    end
    
    fprintf(fid, '\n');
end


fid = fopen('eps_mu_22.txt','w');
for n=1:j-1
    for m=1:i-1
        
    fprintf(fid, '%f', eps_22_matrix(n,m));
    fprintf(fid,'\t');
    end
    
    fprintf(fid, '\n');
end

% virtual_determinant = sqrt((g11_virtual.*g22_virtual - g12_virtual.*g21_virtual));
% physical_determinant = sqrt((g11_physical.*g22_virtual - g12_physical.*g21_physical)); 
