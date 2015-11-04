%Riemann Tensor
%Timothy M. McManus Jr.

clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing Pointwise .wrl file      %
% and extracting desired cooridnates %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load Surface
fid = fopen('gaussian_pulse_complete_201x201.wrl','r'); % Pointwise export file
textscan(fid, '%*s %*[^\n]', 5); % Skip header text
data = textscan(fid, '%f %f %f %*[^\n]', 'delimiter','','MultipleDelimsAsOne', 1, 'CollectOutput', 1); % read numeric data
fclose(fid);
X = data{1}(:,1);
Y = data{1}(:,2);
Z = data{1}(:,3);

i=201; %Points in the i-direction 
j=201; %number of POINTS in the j-direction

%determining the number of vectors given in the i,j and k directions
i_vectors = (i-1).*(j);
j_vectors = (i).*(j-1);
cells=(i-1).*(j-1);

i_k_cells=(i-2).*(j);
j_k_cells=i.*(j-2);
k_cells=(i-2).*(j-2);

%Loading the Points that will be inputted into the FFF:
%X_fff=X; Y_fff=Y; Z_fff=Z;
 
%Removing the Edges from X_fff,Y_fff and Z_fff s.t. the metric tensors can
%be compared directly with those determined by the central difference method
[X_central,Y_central,Z_central]=central_points(X,Y,Z,i,j);

% X_fff(i:i:end)=[]; X_fff(((i-1).*(j-1))+1:((i-1).*j))=[];
% Y_fff(i:i:end)=[]; Y_fff(((i-1).*(j-1))+1:((i-1).*j))=[];
% Z_fff(i:i:end)=[]; Z_fff(((i-1).*(j-1))+1:((i-1).*j))=[];
 
% %Inputting X_fff, Y_fff and Z_fff into the Correct Parametric First
% %Fundamental Forms
% 
% %Hyperbolic Paraboloid
% E = 1+(Y_central).^2; F=(X_central.*Y_central); G = 1+(X_central).^2;
% L=zeros(length(X_central),1);M=sqrt(1+(X_central).^2 + (Y_central).^2); N=zeros(length(X_central),1);
% K=-1./((1+(X_central).^2 + (Y_central).^2)).^2;
% hx=Y_central; hy=X_central; hxx=zeros(length(X_central),1); hyy=hxx; hxy= ones((length(X_central)),1);
% K_orient = (hxx.*hyy - (hxy).^2)./(1+(hx).^2 + (hy).^2).^2;

% %%Prolate Spheroid
% % a = 1; b = 2; c = 3; v = acos(Z_fff./c); u = atan((a.*Y_fff)./(b.*X_fff));
% % E = (b.^2.*((cos(u)).^2) + a.^2.*((sin(u)).^2)).*(sin(u)).^2;
% % F = (b.^2 - a.^2).*(cos(u).*sin(u).*cos(v).*sin(v));
% % G = (a.^2.*(cos(u)).^2 +b.^2.*(sin(u)).^2).*(cos(v)).^2 + c.^2.*(sin(v)).^2;

% %%Shoe Surface
% E = 1 + (X_central).^4; F = -1.*((X_central).^2).*(Y_central); G = 1 + (Y_central).^2;
% L = (2.*X_central)./(sqrt(1+(X_central).^4 + (Y_central).^2));M=zeros(size(X_central),1); N=-1./(sqrt(1+(X_central).^4 + (Y_central).^2));
% K = -(2.*(X_central))./(1+(X_central).^4 + (Y_central).^2).^2;
% hx = (X_central).^2; hy=(-1.*Y_central); hxx=(2.*(X_central)); hyy=-1.*ones(length(X_central),1); hxy=zeros(length(X_central),1);

%%Monkey Saddle
% E = 1 +9.*(((X_central).^2 - (Y_central).^2)).^2; F = -18.*(X_central.*Y_central).*((X_central).^2 - (Y_central).^2); G = 1 + 36.*((X_central.^2).*(Y_central).^2);
% L = (6.*(X_central))./(sqrt(1+9.*((X_central).^2 + (Y_central).^2).^2)); M=(-6.*(Y_central))./(sqrt(1+9.*((X_central).^2 + (Y_central).^2).^2));N = (-6.*(X_central))./(sqrt(1+9.*((X_central).^2 + (Y_central).^2).^2));
% K = (-36.*((X_central).^2 + (Y_central).^2))./(1+9.*((X_central).^2 + (Y_central).^2).^2).^2; 
% K_other = (L.*N - M.^2)./(E.*G - F.^2);

%%Flat Plane
% E = ones(length(X_central),1); F=zeros(length(X_central),1); G = ones(length(X_central),1); L=zeros(length(X_central),1); M=zeros(length(X_central),1); N=zeros(length(X_central),1);
% K=zeros(length(X_central),1);

%Tilted Plane
% E = 2.*ones(length(X_central),1); F=zeros(length(X_central),1); G = ones(length(X_central),1); L=zeros(length(X_central),1); M=zeros(length(X_central),1); N=zeros(length(X_central),1);
% K=zeros(length(X_central),1);hx=ones(length(X_central),1);  hy=zeros(length(X_central),1); hxy=zeros(length(X_central),1); hxx = zeros(length(X_central),1); hyy=zeros(length(X_central),1);

%%PsuedoSphere
% a=1;
% E=(tanh(X_central)).^2; F=zeros(length(X_central),1); G=(sech(X_central)).^2; e = -sech(X_central).*tanh(X_central); f = zeros(length(X_central),1); g = sech(X_central).*tanh(X_central);
% K=-1./a.^2;

%%Sphere
% a = 1;  K=(1./a^2).*ones(length(X_central));  
% for n=1:(length(X_central)-1)
%     h(n)=X_central(n+1)-X_central(n);
% end
% 
% for n=1:(length(Y_central)-(i-2))
%     k(n)=Y_central(n+(i-2))-Y_central(n);
% end

% for n=1:length(Y_central)
%     norm_vector(n,:)=[-1.*hx(n),-1.*hy(n),1];
%     unit_norm_vector(n,:)=norm_vector(n,:)./norm(norm_vector(n,:));
% end
% 
% n_theo = unit_norm_vector;

%functionalized creation a_virtual
[a,a_vector] = a_function(X,Y,Z,i,j);
%functionalized creation of d_virtual
[d,d_vector] = d_function(X,Y,Z,i,j);
%assigning function
[i_dist,i_vector,j_dist,j_vector,c,c_vector,b,b_vector]=assignment(a,a_vector,d,d_vector);
%prooning a,b,c,d to agree with expected result
a_optical = a;
a_optical(1:i-1)=[];
a_optical(1:i-1:end)=[];
a_optical(i-2:i-2:end)=[];
a_optical(end-(i-4):end)=[];
a_optical_mat = (reshape(a_optical,i-3,j-2))';
curved_path_length_horizontal = sum(a_optical_mat((i-1)./2,:));

d_optical = d;
d_optical(1:i)=[];
d_optical(1:i:end)=[];
d_optical(i-1:i-1:end)=[];
d_optical(end-(i-3):end)=[];
d_optical_mat = (reshape(d_optical,j-2,i-3))';
curved_path_length_vertical = sum(d_optical_mat(:,(j-1)./2));


[a,a_vector,c,c_vector,b,b_vector,d,d_vector] = prooning(a,a_vector,c,c_vector,b,b_vector,d,d_vector,i,j);



a_spacing = a_vector; 
a_spacing(1:i-1)=[]; 
a_spacing(1:i:end)=[]; 
a_spacing(i-2:i-2:end)=[];

% length(i_virtual)
% i_virtual((((i-1).*j) - (i-2)):j.*(i-1))=[]; %Removing the top
% i_virtual(i-1:i-1:end)=[];
% length(i_virtual)


% %note that the cross product is going to depend on intial mesh
% orientation
n_vector=cross(a_vector,d_vector);

for n=1:length(n_vector)
    norm_virtual(n)=norm(n_vector(n,:));
    norm_n(n,:)=n_vector(n,:)./(norm_virtual(n));
end

norm_n_row1=norm_n(:,1); norm_n_row2=norm_n(:,2); norm_n_row3=norm_n(:,3);
%Populating 
% norm_n_row1(1:i-1)=[];  norm_n_row1(1:i:end)=[];
% norm_n_row2(1:i-1)=[];  norm_n_row2(1:i:end)=[];
% norm_n_row3(1:i-1)=[];  norm_n_row3(1:i:end)=[];

n_geo=[norm_n_row1,norm_n_row2,norm_n_row3];

g11_central=zeros(1,(i.*j)-1);

[h11_central,h12_central_x,hx_calc,f11_central] = partial_x(Z,X,i,j);
%[h11_central,h11_new,h12_central_x,f11_central]=new_partial_x(Z_virtual_new,X_virtual,i,j);

[h22_central,h12_central_y,f22_central]=partial_y(Z,Y,i,j);
%h12_central_raw = h12_central_x_raw.*h12_central_y_raw;

h12_central=h12_central_x.*h12_central_y;
%Populating covariant and covariant tensors

%Metric tensor area calculator
area = sqrt((h11_central).*(h22_central) - (h12_central).^2);


epsxx_new = (1./area).*(h22_central);
epsxy_new = zeros(1,length(h12_central));
epsyy_new = (1./area).*(h11_central);

epsxx_new_matrix = reshape(epsxx_new,i-2,i-2);
epsxy_new_matrix = reshape(epsxy_new,i-2,i-2);
epsyy_new_matrix = reshape(epsyy_new,i-2,i-2);

epsxx_vertical = epsxx_new_matrix(:,100);
epsxy_vertical = epsxy_new_matrix(:,100);
epsyy_vertical = epsyy_new_matrix(:,100);

nxx_vertical = abs(epsxx_vertical);
nxy_vertical = abs(epsxy_vertical);
nyy_vertical = abs(epsyy_vertical);

nxx_matrix = epsxx_new_matrix;
nxy_matrix = epsxy_new_matrix;
nyy_matrix = epsyy_new_matrix;

flat_inc_length = linspace(min(X_central),max(X_central),i-2);
flat_nxx_optical_length_vertical = trapz(flat_inc_length,nxx_vertical);
flat_nyy_optical_length_vertical = trapz(flat_inc_length,nyy_vertical);
flat_nxy_optical_length_vertical = trapz(flat_inc_length,nxy_vertical);

epsxx_horizontal = epsxx_new_matrix(100,:);
epsxy_horizontal = epsxy_new_matrix(100,:);
epsyy_horizontal = epsyy_new_matrix(100,:);

nxx_horizontal = abs(epsxx_horizontal);
nxy_horizontal = abs(epsxy_horizontal);
nyy_horizontal = abs(epsyy_horizontal);

flat_nxx_optical_length_horizontal = trapz(flat_inc_length,nxx_horizontal);
flat_nyy_optical_length_horizontal = trapz(flat_inc_length,nyy_horizontal);
flat_nxy_optical_length_horizontal = trapz(flat_inc_length,nxy_horizontal);
%Variable assignment for COMSOL input purposes
X=linspace(-5,5,201); Y = linspace(-5,5,201);
[x,y] = meshgrid(X,Y);
r = x.^2 + y.^2;
g11_analytic = 1 + 4.*(x.^2).*exp(-2.*(r));
g12_analytic = 4.*x.*y.*exp(-2.*(r));
g22_analytic = 1 + 4.*(y.^2).*exp(-2.*(r));
g = g11_analytic.*g22_analytic - (g12_analytic).^2;

nxx = g22_analytic./sqrt(g);
nyy = g11_analytic./sqrt(g);
nzz = 1./sqrt(g);

nxx_matrix = reshape(nxx,length(X),length(X));
nyy_matrix = reshape(nyy,length(X),length(X));
nzz_matrix = reshape(nzz,length(X),length(X));

figure
imagesc(nxx_matrix)
figure
imagesc(nyy_matrix)
figure
imagesc(nzz_matrix)

[covariant, contravariant] = tensor_pop(h11_central,h12_central,h22_central,i,j);
 
%Placing variables with know correct solutions to determing where error is
%originating from with regards to solution to Quartered Sphere problem

% Real_E = 1 + (sin(Y_central)).^2;
% Real_F = zeros(1,length(Y_central));
% Real_G = ones(1,length(Y_central));
% 
% h11_central=Real_E;
% h12_central=Real_F;
% h22_central=Real_G;



%X needs to be stripped to agree with the "tightened" mesh
%Calculating Partial Derivatives of metric tensors with respect to x
[pdxh11,pdxh12,pdxh22,X_f1] = pdxh(X,i,h11_central,h12_central,h22_central);

%Y needs to be stripped to agree with the "tightened" mesh
%Calculation Partial Derivatives of metric tensors with respect to y
[pdyh11,pdyh12,pdyh22,Y_f1,pdyh11_raw, pdyh12_raw, pdyh22_raw] = pdyh(Y,j,i,h11_central,h12_central,h22_central);

%Populating Analytic Solutions to Determine Absolute Error
[X_f2,Y_f2,pdxh11_true,pdxh12_true,pdxh22_true,pdyh11_true,pdyh12_true,pdyh22_true,pdxh11_error,pdxh12_error,pdxh22_error,pdyh11_error,pdyh12_error,pdyh22_error] = metric_partial_derivative_error(X_f1,Y_f1,i,pdxh11,pdxh12,pdxh22,pdyh11,pdyh12,pdyh22);

%Populating Chrisotoffel symbols of the first kind (Numeric Solution)
[gamma_111,gamma_112,gamma_221,gamma_222,gamma_121,gamma_122] = Chris_first_kind(pdxh11,pdxh12,pdyh12,pdyh22,pdyh11,pdxh22);

%Populating analytic solutions to determine absolute error in calculation
%of Christoffel symbols of the first kind
[gamma_111_error,gamma_112_error,gamma_221_error,gamma_222_error,gamma_121_error,gamma_122_error,gamma_111_true,gamma_112_true,gamma_221_true,gamma_222_true,gamma_121_true,gamma_122_true] = Chris_first_kind_error(pdxh12_true, pdxh11_true, pdyh12_true, pdxh22_true,pdyh22_true,pdyh11_true,gamma_111,gamma_112,gamma_221,gamma_222,gamma_121,gamma_122);

%Populating Christoffel symbols of the second kind
[gamma_1_11,gamma_1_21,gamma_1_22,gamma_2_11,gamma_2_12,gamma_2_22,gamma_1_11_true,gamma_1_21_true,gamma_1_22_true,gamma_2_11_true,gamma_2_12_true,gamma_2_22_true,contra_h11,contra_h12,contra_h22] = Chris_second_kind(i,j,h11_central,h12_central,h22_central,gamma_111,gamma_112,gamma_221,gamma_222,gamma_121,gamma_122,X_f2,Y_f2,gamma_111_true,gamma_121_true,gamma_221_true,gamma_112_true,gamma_122_true,gamma_222_true);

%Determining absolute error in numeric calculation of the christoffel
%symbol of the second kind
[gamma_1_11_error,gamma_1_21_error,gamma_1_22_error,gamma_2_11_error,gamma_2_12_error,gamma_2_22_error] = Chris_second_kind_error(gamma_1_11,gamma_1_21,gamma_1_22,gamma_2_11,gamma_2_12,gamma_2_22,gamma_1_11_true,gamma_1_21_true,gamma_1_22_true,gamma_2_11_true,gamma_2_12_true,gamma_2_22_true);


%Solving partial derivatives of Christoffel symbols of the second kind
[pdx_gamma_2_21, pdx_gamma_1_21,pdx_gamma_2_22,pdx_gamma_1_22] = pdx_gamma(i,X_f2,gamma_2_12,gamma_1_21,gamma_2_22,gamma_1_22);

[pdy_gamma_2_11,pdy_gamma_1_11,pdy_gamma_2_12,pdy_gamma_1_21] = pdy_gamma(i,Y_f2,gamma_2_11,gamma_1_11,gamma_2_12,gamma_1_21);

%Populating Riemman Tensor of the secod kind
%Strip Christoffel symbols of the second kind such that their length =
%length of pdx_gammas

[R_2_112,R_1_112,R_2_221,R_1_212,R_2_121,R_1_121,R_2_212,R_1_221] = Riemann_2nd_kind(gamma_2_12,gamma_1_21,gamma_2_22,gamma_1_22,gamma_2_11,gamma_1_11,i,j,pdx_gamma_2_21,pdy_gamma_2_11,pdx_gamma_1_21,pdy_gamma_1_11,pdy_gamma_2_12,pdx_gamma_2_22,pdx_gamma_1_22,pdy_gamma_1_21);

%Populating Riemann Tensor of the second kind 
[R_2_112_true,R_1_112_true,R_2_221_true,R_1_212_true,R_2_121_true,R_1_121_true,R_2_212_true,R_1_221_true,X_f3,Y_f3]=Riemann_2nd_kind_true(X_f2,Y_f2,gamma_1_21_true,gamma_2_11_true,gamma_2_12_true,gamma_1_22_true,i,j);

%Evaluating Absolute Error of Riemann_2nd_kind_Error
[R_2_112_error, R_1_112_error,R_2_221_error,R_1_212_error] = Riemann_2nd_kind_error(R_2_112_true,R_1_112_true,R_2_221_true,R_1_212_true,R_2_112,R_1_112,R_2_221,R_1_212);


%Calculating Ricci tensor of the second kind and Ricci Scalar
[R_11,R_12,R_21,R_22,R] = Ricci_2nd_kind(R_2_121,R_1_112,R_1_212,contra_h11,contra_h12,contra_h22,i);

%Calculating True Ricci tensor of the second kind and the absolute error
[R_true,R_error] = Ricci_2nd_kind_error(R_2_121_true, R_1_112_true,R_1_212_true,X_f3,Y_f3,R);

%Comparing with Gaussian Curvature assuming Monge Patch
[K,K_R_error_true,K_R_error] = Curvature_comparator(R,R_true,X_f3,Y_f3);

K_matrix = reshape(K,(i-6),(i-6));
R_matrix = reshape(R,(i-6),(i-6));

[f12_central] = xy_partial(Z,X,Y,i,j);

K_orient_calc=(f11_central.*f22_central - (f12_central).^2)./(1+(h12_central_x).^2 + (h12_central_y).^2).^2;

% K_no_error=((f11_central.*f22_central)-(f12_central).^2)./((h12_central_x).*(h12_central_y) - ((h11_central).*(h22_central))).^2;


second_FFF=zeros((i-2).*(j-2),3);
f11_calc=second_FFF; f12_calc=second_FFF; f22_calc=second_FFF;

for n=1:(i-2).*(j-2)
    f11_calc(n,3)=f11_central(n);
    f12_calc(n,3)=f12_central(n);
    f22_calc(n,3)=f22_central(n);
end

%Calculating Elemental Area based on parallelogram area
% elemental_area = sqrt(h11_central.*h22_central - (h12_central).^2);
% 
% %Calculating SFF coefficients using normal vectors and SPD
% [L_normal,M_normal,N_normal,K_normal]=geometric_SFF(f11_central,f22_central,f12_central,h12_central_x,h12_central_y,norm_n,h11_central,h22_central,h12_central);
% 
% %Calculating SFF coefficeints using proposed formalism for Gaussian
% %Curvature (K)
% [L_det,M_det,N_det,K_det]=algebraic_SFF(f11_central,h11_central,h22_central,h12_central,f12_central,f22_central);
% E_calc = h11_central;
% F_calc = h12_central;
% G_calc = h22_central;


%  for n=1:length(f11_central)
%      E_absolute_error(n)=E_calc(n)-E(n);
%      F_absolute_error(n)=F_calc(n)-F(n);
%      G_absolute_error(n)=G_calc(n)-G(n);
%      L_absolute_error(n)=L_det(n)-L(n);   
%      M_absolute_error(n)=M_det(n)-M(n);
%      N_absolute_error(n)=N_det(n)-N(n);
%  end




% figure
% subplot(2,1,1)
% x=1:length(h11_central);
% plot(x,E,'r',x,h11_central,'b-')
% xlabel('index','FontSize',11,'Fontweight','bold')
% xlim([1 length(h11_central)])
% legend('E','E_c_a_l_c')
% grid on
% title(sprintf('E and E_c_a_l_c (%dx%d Atlas)',i-3,j-3),'FontSize',13,'Fontweight','bold')
% subplot(2,1,2)
% plot(E_absolute_error)
% xlabel('index','FontSize',11,'Fontweight','bold')
% xlim([1 length(h11_central)])
% title(sprintf('E Absolute Error (%dx%d Atlas)',i-3,j-3),'FontSize',13,'Fontweight','bold')
% grid on

% 
% figure
% subplot(2,1,1)
% x=1:length(f11_central);
% plot(x,L,'r',x,f11_central,'b-')
% xlabel('index','FontSize',11,'Fontweight','bold')
% xlim([1 length(h11_central)])
% legend('L','L_c_a_l_c')
% grid on
% title(sprintf('L and L_c_a_l_c (%dx%d Atlas)',i-3,j-3),'FontSize',13,'Fontweight','bold')
% subplot(2,1,2)
% plot(L_absolute_error)
% xlabel('index','FontSize',11,'Fontweight','bold')
% xlim([1 length(h11_central)])
% title(sprintf('L Absolute Error (%dx%d Atlas)',i-3,j-3),'FontSize',13,'Fontweight','bold')
% grid on

% figure
% plot(E_absolute_error,'b','Linewidth',2)
% title(sprintf('E Absolute Error (%dx%d Atlas)',i-3,j-3),'FontSize',13,'Fontweight','bold')
% xlabel('Patch index','FontSize',11,'Fontweight','bold')
% xlim([1 length(E_absolute_error)])
% 



for n=1:k_cells
    E_prime(n) = (1./((h11_central(n).*h22_central(n) - (h12_central(n)).^2))).*h22_central(n);
    F_prime(n)= (1./((h11_central(n).*h22_central(n) - (h12_central(n)).^2))).*-.1*h12_central(n);
    G_prime(n)= (1./((h11_central(n).*h22_central(n) - (h12_central(n)).^2))).*h11_central(n);
end



% E_prime_new_matrix =reshape(h11_new,[],j-1);
% F_prime_new_matrix=reshape(h12_new,[],j-1);
% G_prime_new_matrix=reshape(h22_new,[],j-1);

a_spacing_matrix = reshape(a_spacing,[],j-2);
% K_matrix=reshape(K,[],j-2);
% K_no_error_matrix=reshape(K_no_error,[],j-2);



% K_calc_matrix=reshape(K_det,[],j-2);
% figure
% imagesc(flipud(K_calc_matrix'))
% title(sprintf('K_c_a_l_c (%dx%d Surface Mesh)',i-2,j-2),'FontSize',13,'Fontweight','bold')
% xlim([1 i-2])
% ylim([1 j-2])
% xlabel('X(mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-2,3),'XTickLabel',(linspace(min(X_central),max(X_central),3)),'Fontweight','bold')
% ylabel('Y (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-2,3),'YTickLabel',(linspace(max(Y_central),min(Y_central),3)))
% colorbar
% 
% K_error = K-K_det';
% 
% x=1:length(K_error);
% figure
% plot(x,K_error,'b','Linewidth',2)
% title(sprintf('Absolute Error(%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlabel('Point index','FontSize',11,'Fontweight','bold')
% %ylabel('degrees','FontSize',11,'Fontweight','bold')
% xlim([1 length(K_error)])


% figure
% imagesc(flipud(E_prime_new_matrix'))
% title(sprintf('E (%dx%d Surface Mesh)',i-1,j-1),'FontSize',13,'Fontweight','bold')
% xlim([1 i-1])
% ylim([1 j-1])
% xlabel('Pos. along x-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'XTick',linspace(1,j-1,11),'XTickLabel',(linspace(min(X_virtual),max(X_virtual),11)),'Fontweight','bold')
% ylabel('Pos. along y-axis (mm)','FontSize',11,'Fontweight','bold')
% set(gca,'YTick',linspace(1,i-1,11),'YTickLabel',(linspace(max(Y_virtual),min(Y_virtual),11)))
% colorbar

