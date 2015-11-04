%TEC_3D_2p0_COMSOL_centroid
%Timothy M. McManus Jr.
%Last edited: 26-8-2014
clear all
close all
clc

tic
%MATLAB Live Link
import com.comsol.model.*
import com.comsol.model.util.*

%Creating a hexadron
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing Pointwise .uns file      %
% and extracting desired cooridnates %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Virtual Space
fid = fopen('x_y_virtual_block.uns','r'); % Pointwise export file
textscan(fid, '%*s %*[^\n]', 16); % Skip header text
Virtual_data = textscan(fid, '%f %f %f %*[^\n]', 'delimiter','','MultipleDelimsAsOne', 1, 'CollectOutput', 1); % read numeric data
fclose(fid);

%NOTE THAT Y and Z Are swapped for both spaces

X_virtual = Virtual_data{1}(:,1);
Z_virtual = Virtual_data{1}(:,2);
Y_virtual = Virtual_data{1}(:,3);
Virtual = [X_virtual,Y_virtual,Z_virtual];

%Load Physical Space
fid = fopen('x_y_physical_block.uns','r');
textscan(fid,'%*s %*[^\n]', 16);
Physical_data = textscan(fid, '%f %f %f %*[^\n]', 'delimiter','','MultipleDelimsAsOne', 1, 'CollectOutput', 1);
fclose(fid);

%NOTE THAT Y and Z Are swapped for both spaces
X_physical = Physical_data{1}(:,1);
Z_physical = Physical_data{1}(:,2);
Y_physical = Physical_data{1}(:,3);
Physical = [X_physical,Y_physical,Z_physical];
t_import = toc;
disp(sprintf('t_import: %d sec', t_import));

%ensure that Z_physical and Z_virtual are being imported in the same
%direction, i.e. (z1 to z2)

tic;
i=101; %number of POINTS in the horizontal direction(i-direction) (major axis)
j=101;%number of POINTS in the depth direction (j-direction)
k=101;%number of POINTS in the vertical direction(k-direction)
eps_virtual = 2.6;

%determining the number of vectors given in the i,j and k directions
i_vectors = (i-1)*(j)*(k);
j_vectors = (i)*(j-1)*(k);
k_vectors = (i)*(j)*(k-1);
cuboids = (i-1)*(j-1)*(k-1);

parfor n = 1:(i*j*k)-1
    a_physical(n)=sqrt((X_physical(n)-X_physical(n+1)).^2 + (Y_physical(n)-Y_physical(n+1)).^2 + (Z_physical(n) - Z_physical(n+1)).^2);
    a_virtual(n)= sqrt((X_virtual(n)-X_virtual(n+1)).^2 + (Y_virtual(n)-Y_virtual(n+1)).^2 + (Z_virtual(n) - Z_virtual(n+1)).^2);
    a_physical_vector(n,:)=[X_physical(n+1)-X_physical(n),Y_physical(n+1)-Y_physical(n), +Z_physical(n+1)-Z_physical(n)];
    a_virtual_vector(n,:)=[X_virtual(n+1)-X_virtual(n),Y_virtual(n+1)-Y_virtual(n), Z_virtual(n+1)-Z_virtual(n)];
    a_physical_n = n;
end

parfor n=1:(i*j*k)-i
    c_physical(n)=sqrt((X_physical(n) - X_physical(n+i)).^2 +(Y_physical(n)-Y_physical(n+i)).^2+(Z_physical(n)-Z_physical(n+i)).^2);
    c_virtual(n)=sqrt((X_virtual(n) - X_virtual(n+i)).^2 +(Y_virtual(n)-Y_virtual(n+i)).^2+(Z_virtual(n)-Z_virtual(n+i)).^2);
    c_physical_vector(n,:)=[X_physical(n+i)-X_physical(n),Y_physical(n+i)-Y_physical(n), +Z_physical(n+i)-Z_physical(n)];
    c_virtual_vector(n,:)=[X_virtual(n+i)-X_virtual(n),Y_virtual(n+i)-Y_virtual(n), Z_virtual(n+i)-Z_virtual(n)];
    c_physical_n = n;
end

parfor n=1:(i.*j*k)-(i*j)
    b_physical(n)=sqrt((X_physical(n)-X_physical(n+(i.*j))).^2 +(Y_physical(n)- Y_physical(n+(i.*j))).^2 + (Z_physical(n)-Z_physical(n+(i.*j))).^2);
    b_virtual(n)=sqrt((X_virtual(n)-X_virtual(n+(i.*j))).^2 +(Y_virtual(n)- Y_virtual(n+(i.*j))).^2 + (Z_virtual(n)-Z_virtual(n+(i.*j))).^2);
    b_physical_vector(n,:)=[X_physical(n+(i.*j))-X_physical(n),Y_physical(n+(i.*j))-Y_physical(n), +Z_physical(n+(i.*j))-Z_physical(n)];
    b_virtual_vector(n,:)=[X_virtual(n+(i.*j))-X_virtual(n),Y_virtual(n+(i.*j))-Y_virtual(n), Z_virtual(n+(i.*j))-Z_virtual(n)];
    b_physical_n = n;
end

%Removing artefacts that occur as a result of counting technique
a_physical(i:i:end)=[];
a_virtual(i:i:end)=[];
a_physical_vector(i:i:end,:)=[];
a_virtual_vector(i:i:end,:)=[];
%Only valid for when i=j=k;

remove=[(j*(i-1))+1:(i*j)];
parfor n=1:k-2
    remove_added(n,:)=(n.*(i.*j))+remove;
end
remove=[remove;remove_added];
c_physical(remove)=[];
c_virtual(remove)=[];
c_physical_vector(remove,:)=[];
c_virtual_vector(remove,:)=[];
%Only valid for when i=j=k;

a1_physical = a_physical;
a1_virtual =a_virtual;
a1_physical_vector = a_physical_vector;
a1_virtual_vector = a_virtual_vector;
a1_seed = (((i-1)*(j-1)+1):(i-1)*j);

parfor n=1:k-2
    a1_remove(n,:)=n.*((j.*(i-1)))+a1_seed;
end

a1_remove = [a1_seed;a1_remove];
a1_back = (((i-1).*j)*(k-1))+1:(k.*(i-1)*j);
a1_back=reshape(a1_back,[],(i-1));
a1_remove =[a1_remove;a1_back];
a1_remove = reshape(a1_remove,[],1);
a1_physical = a1_physical(setdiff(1:length(a1_physical),a1_remove));
a1_virtual = a1_virtual(setdiff(1:length(a1_virtual),a1_remove));

a1_remove=flipud(sort(a1_remove));
for n=1:length(a1_remove)
    a1_physical_vector(a1_remove(n),:)=[];
    a1_virtual_vector(a1_remove(n),:)=[];
end

a2_physical = a_physical;
a2_virtual =a_virtual;
a2_seed = (1:(i-1));
a2_physical_vector = a_physical_vector;
a2_virtual_vector = a_virtual_vector;

parfor n=1:k-2
    a2_remove(n,:)=n.*((j.*(i-1)))+a2_seed;
end

a2_remove = [a2_seed;a2_remove];
a2_back = (((i-1)*j)*(k-1))+1:(k.*(i-1)*j);
a2_back=reshape(a2_back,[],(i-1));
a2_remove =[a2_remove;a1_back];
a2_remove = reshape(a2_remove,[],1);
a2_physical = a2_physical(setdiff(1:length(a2_physical),a2_remove));
a2_virtual = a2_virtual(setdiff(1:length(a2_virtual),a2_remove));

a2_remove=flipud(sort(a2_remove));
for n=1:length(a2_remove)
    a2_physical_vector(a2_remove(n),:)=[];
    a2_virtual_vector(a2_remove(n),:)=[];
end

a3_physical = a_physical;
a3_virtual =a_virtual;
a3_physical_vector = a_physical_vector;
a3_virtual_vector = a_virtual_vector;

a3_physical(1:(i-1)*j)=[];
a3_first = [1:(i-1)*j];
a3_first = reshape(a3_first,[],i-1);
a3_virtual(1:(i-1)*j)=[];
a3_seed = ((i-1)*j) - (i-2) : (i-1)*j;
parfor n = 1:k-1
    a3_remove(n,:)=n.*a3_seed;
end

a3_remove_total = [a3_remove;a3_first];
a3_remove_total = reshape(a3_remove_total,[],1);
a3_remove_total = flipud(sort(a3_remove_total));

for n=1:length(a3_remove_total)
    a3_physical_vector(a3_remove_total(n),:)=[];
    a3_virtual_vector(a3_remove_total(n),:)=[];
end

a3_remove = reshape(a3_remove,[],1);
a3_physical = a3_physical(setdiff(1:length(a3_physical),a3_remove));
a3_virtual = a3_virtual(setdiff(1:length(a3_virtual),a3_remove));

%a4, by definition has no place on the front and bottom face of the
%structure

a4_physical = a_physical;
a4_virtual =a_virtual;
a4_physical_vector = a_physical_vector;
a4_virtual_vector = a_virtual_vector;
%a4_seed = (1:(i-1));

%Removing the front face
a4_physical(1:(i-1)*j)=[];
a4_virtual(1:(i-1)*j)=[];
a4_first = [1:(i-1)*j];

%Removing the bottom face
a4_seed = 1:i-1;
parfor n=1:k-1
    a4_remove(n,:)=(a4_seed)+(n-1)*((i-1)*j);
end

a4_first = reshape(a4_first,[],i-1);
a4_remove_total = [a4_remove;a4_first];
a4_remove_total = reshape(a4_remove_total,[],1);
a4_remove_total = flipud(sort(a4_remove_total));

for n=1:length(a4_remove_total)
    a4_physical_vector(a4_remove_total(n),:)=[];
    a4_virtual_vector(a4_remove_total(n),:)=[];
end

a4_remove = reshape(a4_remove,[],1);
a4_physical = a4_physical(setdiff(1:length(a4_physical),a4_remove));
a4_virtual = a4_virtual(setdiff(1:length(a4_virtual),a4_remove));

% To this point a-edge lengths and a-oriented vectors appear to be correct

b1_physical = b_physical;
b1_virtual = b_virtual;
b1_physical_vector = b_physical_vector;
b1_virtual_vector = b_virtual_vector;

%Remove right edges, of which b1 is not defined along the "right face"

b1_right_seed = [i:i:i*j];
parfor n = 1:(k-1)
    b1_right(n,:)=(n-1)*(i*j) + b1_right_seed;
end
b1_right = reshape(b1_right,[],1);

b1_top_seed = [((i-1)*j)+1:(i*j)-1];
parfor n = 1:(k-1)
    b1_top(n,:)=(n-1)*(i*j)+b1_top_seed;
end
b1_top = reshape(b1_top,[],1);
b1_remove=sort([b1_right;b1_top]);

b1_physical = b1_physical(setdiff(1:length(b1_physical),b1_remove));
b1_virtual = b1_virtual(setdiff(1:length(b1_virtual),b1_remove));
b1_remove = flipud(b1_remove);
for n=1:length(b1_remove)
    b1_physical_vector(b1_remove(n),:)=[];
    b1_virtual_vector(b1_remove(n),:)=[];
end

%b1 edges and vectores appear to be behaving properly now

b2_physical = b_physical;
b2_virtual = b_virtual;
b2_physical_vector = b_physical_vector;
b2_virtual_vector = b_virtual_vector;

b2_right_seed = [i:i:i*(j-1)];
parfor n = 1:i-1
    b2_right(n,:) = (n-1)*(i*j)+b2_right_seed;
end
b2_right = reshape(b2_right,[],1);
b2_top_seed = [(i*(j-1))+1:i*j];
 parfor n = 1:i-1
     b2_top(n,:) = (n-1)*(i*j)+b2_top_seed;
 end
 b2_top = reshape(b2_top,[],1);
 b2_remove=sort([b2_top;b2_right]);
 b2_physical = b2_physical(setdiff(1:length(b2_physical),b2_remove));
 b2_virtual = b2_virtual(setdiff(1:length(b2_virtual),b2_remove));

 b2_remove = flipud(b2_remove);

 for n=1:length(b2_remove)
    b2_physical_vector(b2_remove(n),:)=[];
    b2_virtual_vector(b2_remove(n),:)=[];
 end

%b2 edges and vectors appear to be behaving properly

b3_physical = b_physical;
b3_virtual = b_virtual;

b3_physical_vector = b_physical_vector;
b3_virtual_vector = b_virtual_vector;

%Need to remove the right side layer from the structure (while not touching
%the very bottom layer)

parfor n=1:j-1
    b3_seed(1,n)=(n+1).*i;
end

parfor n=1:k-2
    b3_remove(n,:)=n.*(j.*i)+b3_seed;
end

b3_remove=[b3_seed;b3_remove];

bottom_seed = [1:i];
parfor n=1:k-2
    bottom_remove(n,:)= n.*(j.*i)+bottom_seed;
end

bottom_remove = [bottom_seed;bottom_remove];
bottom_remove=reshape(bottom_remove,[],1);
b3_remove=reshape(b3_remove,[],1);
b3_remove_total =[bottom_remove;b3_remove];

b3_remove_total=flipud(sort(b3_remove_total));
for n=1:length(b3_remove_total)
    b3_physical_vector(b3_remove_total(n),:)=[];
    b3_virtual_vector(b3_remove_total(n),:)=[];
end

b3_physical = b3_physical(setdiff(1:length(b3_physical),b3_remove_total));
b3_virtual = b3_virtual(setdiff(1:length(b3_virtual),b3_remove_total));

%b3 edges and vectors appear to be behaving properly

b4_physical = b_physical;
b4_virtual = b_virtual;
b4_physical_vector = b_physical_vector;
b4_virtual_vector = b_virtual_vector;

%This removes all of the edges on the far-left of the structure (except the
%bottom layer)

b4_left_seed = [i+1:i:(i*j)-(i-1)];

parfor n = 1:i-1
    b4_left(n,:)=(n-1)*(i*j)+b4_left_seed;
end
b4_left = reshape(b4_left,[],1);

b4_bottom_seed =[1:i];
parfor n=1:i-1
    b4_bottom(n,:)= (n-1).*((i*j))+b4_bottom_seed;
end
b4_bottom=reshape(b4_bottom,[],1);

b4_remove = sort([b4_bottom;b4_left]);
b4_physical = b4_physical(setdiff(1:length(b4_physical),b4_remove));
b4_virtual = b4_virtual(setdiff(1:length(b4_virtual),b4_remove));

b4_remove = flipud(b4_remove);
for n=1:length(b4_remove)
    b4_physical_vector(b4_remove(n),:)=[];
    b4_virtual_vector(b4_remove(n),:)=[];
end

%b4 edges and vectors appear to be behaving properly


%Defining C edges
c1_physical = c_physical;
c1_virtual = c_virtual;
c1_physical_vector = c_physical_vector;
c1_virtual_vector = c_virtual_vector;

parfor n=1:j-1
    seed(1,n)=n.*i;
end

parfor n=1:k-2
    seed_c1(n,:)=n.*(i.*(j-1))+seed;
end

c1_remove=[seed;seed_c1];
c1_back = i.*(j-1).*(k-1)+1:i.*(j-1).*k;
c1_back=reshape(c1_back,[],(j-1));
c1_remove=[c1_remove;c1_back];
c1_remove = reshape(c1_remove,[],1);
c1_remove = sort(c1_remove);
c1_physical = c1_physical(setdiff(1:length(c1_physical),c1_remove));
c1_virtual = c1_virtual(setdiff(1:length(c1_virtual),c1_remove));

c1_remove=flipud(c1_remove);
for n=1:length(c1_remove)
    c1_physical_vector(c1_remove(n),:)=[];
    c1_virtual_vector(c1_remove(n),:)=[];
end


%c1 edges and vectors appear to be well behaved

c2_physical = c_physical;
c2_virtual = c_virtual;
c2_physical_vector = c_physical_vector;
c2_virtual_vector = c_virtual_vector;

c2_left_seed = [1:i:i*(j-1)];
parfor n = 1:i-1
    c2_left(n,:) = (n-1)*(i-1)*(j) + c2_left_seed;
end
c2_left = reshape(c2_left,[],1);

c2_back_seed = [k*(j-1)*i - (i-1):k*(j-1)*i];
parfor n=1:i-1
    c2_back(n,:)= c2_back_seed - (n-1)*i;
end
c2_back = reshape(c2_back,[],1);

c2_remove = [c2_back;c2_left];
c2_remove = sort(c2_remove);
c2_physical = c2_physical(setdiff(1:length(c2_physical),c2_remove));
c2_virtual = c2_virtual(setdiff(1:length(c2_virtual),c2_remove));

c2_remove=flipud(c2_remove);
for n=1:length(c2_remove)
    c2_physical_vector(c2_remove(n),:)=[];
    c2_virtual_vector(c2_remove(n),:)=[];
end

%c2 edges and vectors appaer to be behaving properly

c3_physical = c_physical;
c3_virtual = c_virtual;
c3_physical_vector = c_physical_vector;
c3_virtual_vector = c_virtual_vector;

%using linspace for the indice is lazy, and dangerous
seed = linspace(i.*j,i.*(j-2)+i.*j,j-1);

%This is an ugly seed - think about stripping it down
c3_left_seed = [(i*(j-1))+1:i:(i-1)*(j-1)+((i-1)*j)];
parfor n=1:k-1
    c3_left(n,:)=(n-1)*((i-1)*j)+c3_left_seed;
end
c3_left = reshape(c3_left,[],1);
c3_front_seed = [1:1:i];
parfor n = 1:j-1
    c3_front(n,:) = (n-1)*(i)+c3_front_seed;
end
c3_front = reshape(c3_front,[],1);
c3_remove = sort([c3_front;c3_left]);
c3_physical = c3_physical(setdiff(1:length(c3_physical),c3_remove));
c3_virtual = c3_virtual(setdiff(1:length(c3_virtual),c3_remove));

c3_remove=flipud(c3_remove);
for n=1:length(c3_remove)
    c3_physical_vector(c3_remove(n),:)=[];
    c3_virtual_vector(c3_remove(n),:)=[];
end

%c3 edges and vectors appear to be behaving properly
c4_physical = c_physical;
c4_virtual = c_virtual;
c4_physical_vector = c_physical_vector;
c4_virtual_vector = c_virtual_vector;

c4_left_seed =  [(i*(j-1))+1:i:2*(i*(j-1))];
    
parfor n=1:k-1
    c4_left(n,:) = (n-1)*(i*(j-1)) + c4_left_seed;
end
c4_left = reshape(c4_left,[],1);
c4_front_seed = [1:i];

parfor n = 1:k-1
    c4_front(n,:) = (n-1)*(i) + c4_front_seed;
end
c4_front = reshape(c4_front,[],1);
c4_remove = sort([c4_left;c4_front]);
c4_physical = c4_physical(setdiff(1:length(c4_physical),c4_remove));
c4_virtual = c4_virtual(setdiff(1:length(c4_virtual),c4_remove));

c4_remove=flipud(c4_remove);
for n=1:length(c4_remove)
    c4_physical_vector(c4_remove(n),:)=[];
    c4_virtual_vector(c4_remove(n),:)=[];
end

%c4 edges and vectors appear to be behaving properly


%Hexahedron Volume Calcuator

% for n=1:(i.*j.*k)-(i.*(j+1)+1)
%     diag1_physical_vector(n,:)=[X_physical(n+(i.*(j+1)+1))-X_physical(n),Y_physical(i.*(j+1)+1+n)-Y_physical(n),Z_physical(i.*(j+1)+1+n)-Z_physical(n)];
%     diag1_virtual_vector(n,:)=[X_virtual(i.*(j+1)+1+n)-X_virtual(n),Y_virtual(i.*(j+1)+1+n)-Y_virtual(n),Z_virtual(i.*(j+1)+1+n)-Z_virtual(n)];
% end
%
%
% right_seed=linspace(i,i.*(j-2),j-2);
% right_remove=zeros(k-2,length(right_seed));
% for n=1:k-2
%     right_remove(n,:)= n.*(j.*i)+right_seed;
% end
%
% remade_right_remove=[right_seed;right_remove];
% remade_right_remove=reshape(remade_right_remove,[],1);
% top_seed=[j.*(i-1)-1:(j.*i)];
% top_remove=zeros(k-3,length(top_seed));
%
% for n=1:k-3
%     top_remove(n,:)= n.*((j.*i))+top_seed;
% end
% remade_top_remove=[top_seed;top_remove];
% remade_top_remove=reshape(remade_top_remove,[],1);
%
% remove_total=[remade_right_remove;remade_top_remove];
% remove_total=flipud(sort(remove_total));
%
% for n=1:length(remove_total)
%     diag1_physical_vector(remove_total(n),:)=[];
%     diag1_virtual_vector(remove_total(n),:)=[];
% end

r_16_physical_vector = b1_physical_vector+c3_physical_vector;
r_16_virtual_vector = b1_virtual_vector+c3_virtual_vector;
r_25_physical_vector =c3_physical_vector+ (-1.*(b3_physical_vector));
r_25_virtual_vector = c3_virtual_vector+(-1.*(b3_virtual_vector));

S_1562_physical = (0.5).*(cross(r_16_physical_vector,r_25_physical_vector));
S_1562_virtual = (0.5).*(cross(r_16_virtual_vector,r_25_virtual_vector));

r_18_physical_vector = a1_physical_vector+c2_physical_vector;
r_18_virtual_vector = a1_virtual_vector+c2_virtual_vector;
r_54_physical_vector = (-1.*(c1_physical_vector))+a1_physical_vector;
r_54_virtual_vector = (-1.*(c1_virtual_vector))+a1_virtual_vector;

S_1485_physical=(0.5).*(cross(r_18_physical_vector,r_54_physical_vector));
S_1485_virtual=(0.5).*(cross(r_18_virtual_vector,r_54_virtual_vector));

r_13_physical_vector = a1_physical_vector+b2_physical_vector;
r_13_virtual_vector = a1_virtual_vector+b2_virtual_vector;
r_42_physical_vector=-1.*(a1_physical_vector)+b1_physical_vector;
r_42_virtual_vector=-1.*(a1_virtual_vector)+b1_virtual_vector;

S_1234_physical=(0.5).*(cross(r_13_physical_vector,r_42_physical_vector));
S_1234_virtual =(0.5).*(cross(r_13_virtual_vector,r_42_virtual_vector));

diag1_physical_vector=a1_physical_vector+b2_physical_vector+c4_physical_vector;
diag1_virtual_vector=a1_virtual_vector+b2_virtual_vector+c4_virtual_vector;

parfor n=1:(i-1).*(j-1).*(k-1)
hex_volume_physical(n)= (1./3).*dot((S_1234_physical(n,:)+S_1485_physical(n,:)+S_1562_physical(n,:)),diag1_physical_vector(n,:));
hex_volume_virtual(n)= (1./3).*dot((S_1234_virtual(n,:)+S_1485_virtual(n,:)+S_1562_virtual(n,:)),diag1_virtual_vector(n,:));
end

t_calc = toc;
disp(sprintf('t_calc: %d sec', t_calc));

% physical_cell_info = cell(j-1,i-1,k-1);
% virtual_cell_info = cell(j-1,i-1,k-1);
% 
% 
% for m=1:i-1
%     for n=1:j-1
%         for o=1:k-1
%             physical_cell_info{n,m,o}=zeros(3,4);
%         end
%     end
% end

a_physical_form=[a1_physical',a2_physical',a3_physical',a4_physical'];
a_virtual_form=[a1_virtual',a2_virtual',a3_virtual',a4_virtual'];
b_physical_form=[b1_physical',b2_physical',b3_physical',b4_physical'];
b_virtual_form=[b1_virtual',b2_virtual',b3_virtual',b4_virtual'];
c_physical_form=[c1_physical',c2_physical',c3_physical',c4_physical'];
c_virtual_form=[c1_virtual',c2_virtual',c3_physical',c4_virtual'];

%Formatting for cell entry


%Removing artefacts that arise as a consequence of counting algorithm

cuboid_volume1_physical = (a1_physical.*b1_physical.*c1_physical)';
cuboid_volume2_physical = (a1_physical.*b2_physical.*c2_physical)';
cuboid_volume3_physical = (a4_physical.*b4_physical.*c4_physical)';
cuboid_volume4_physical = (a4_physical.*b3_physical.*c3_physical)';

cuboid_volume1_virtual = (a1_virtual.*b1_virtual.*c1_virtual)';
cuboid_volume2_virtual = (a1_virtual.*b2_virtual.*c2_virtual)';
cuboid_volume3_virtual = (a4_virtual.*b4_virtual.*c4_virtual)';
cuboid_volume4_virtual = (a4_virtual.*b3_virtual.*c3_virtual)';

%Preparing points to be exported to COMSOL
physical_hex_info = cell((i*j)*(k-1) - i -1,1);
virtual_hex_info = cell((i*j)*(k-1) - i -1,1);

%Initializing Cells
for n=1:(i*j)*(k-1) - i -1
       physical_hex_info{n,1}=zeros(8,3);
       virtual_hex_info{n,1}=zeros(8,3);
end

for n = 1:(i*j)*(k-1) - i -1
    physical_hex_info{n,1}=[X_physical(n),Y_physical(n),Z_physical(n);X_physical(n+1),Y_physical(n+1),Z_physical(n+1);
        X_physical(n+i),Y_physical(n+i),Z_physical(n+i);X_physical(n+i+1),Y_physical(n+i+1),Z_physical(n+i+1);X_physical(n+i*j),Y_physical(n+i*j),Z_physical(n+i*j);
        X_physical(n+i*j+1),Y_physical(n+i*j+1),Z_physical(n+i*j+1);X_physical(n+i*j+i),Y_physical(n+i*j+i),Z_physical(n+i*j+i);
        X_physical(n+i*j+i+1), Y_physical(n+i*j+i+1), Z_physical(n+i*j+i+1)];
    
    virtual_hex_info{n,1}=[X_virtual(n),Y_virtual(n),Z_virtual(n);X_virtual(n+1),Y_virtual(n+1),Z_virtual(n+1);
        X_virtual(n+i),Y_virtual(n+i),Z_virtual(n+i);X_virtual(n+i+1),Y_virtual(n+i+1),Z_virtual(n+i+1);X_virtual(n+i*j),Y_virtual(n+i*j),Z_virtual(n+i*j);
        X_virtual(n+i*j+1),Y_virtual(n+i*j+1),Z_virtual(n+i*j+1);X_virtual(n+i*j+i),Y_virtual(n+i*j+i),Z_virtual(n+i*j+i);
        X_virtual(n+i*j+i+1), Y_virtual(n+i*j+i+1), Z_virtual(n+i*j+i+1)];
end

%Pruning method
%This first cut appears to work fine (Validated via COMSOL
%mphgeom(model,'geom1')
physical_hex_info(i:i:end)=[];
virtual_hex_info(i:i:end)=[];
top_seed = [(i-1)*(j-1)+1:i*(j-1)];

%DOUBLE CHECK TO MAKE SURE THIS IS VALID
  for n = 1:k-2
      physical_hex_info((n-1)*((i-1)*(j-1))+top_seed)=[];
      virtual_hex_info((n-1)*((i-1)*(j-1))+top_seed)=[];
  end
  
%Determining Centroids of hexahedral cells
  
for n = 1:length(physical_hex_info)
    physical_hex_centroid(n,:)=[sum(physical_hex_info{n}(:,1))/8,sum(physical_hex_info{n}(:,3))/8,sum(physical_hex_info{n}(:,2))/8];
    virtual_hex_centroid(n,:)=[sum(virtual_hex_info{n}(:,1))/8,sum(virtual_hex_info{n}(:,3))/8,sum(virtual_hex_info{n}(:,2))/8];
    hex_perm(n,1)=eps_virtual*(hex_volume_virtual(n)./hex_volume_physical(n));
end
  
%  phys_length=size(physical_hex_info);
%  phys_length=phys_length(1);
%  physical_hex_info = physical_hex_info(setdiff(1:phys_length,top));
 %c4_virtual = c4_virtual(setdiff(1:length(physical_hex_info),c4_remove));
 
% check = size(physical_hex_info);
% check = check(1);
% if check ~= (i-1)*(j-1)*(k-1)
%     error('Error in Cell Array Pruning');
% end 


%IMPORTANT
%For a hexahedron approximately aligned to the coordinate planes, the points in p
%are ordered as follows. The first four points and the last four points projected down
%to the (x, y)-plane defines two negatively oriented quadrangles. The corresponding
%plane for the second quadrangle must lie above the plane of the first quadrant in the
%z direction. 
plot(hex_perm);
interp_function = [physical_hex_centroid,hex_perm];
dlmwrite('x_y_perm_map.txt',interp_function,'delimiter','\t','precision',4);

model = ModelUtil.create('Model2');
model.param.set('D','63.75 [mm]', 'Lens Diameter')
model.param.set('f','10.0 [GHz]', 'central frequency operating');
%Change the length unit (the default is m)
%model.geom('geom1').lengthUnit('mm');
model.geom().create('geom1',3);
hex_size = size(physical_hex_info);
hex_size = hex_size(1);

for  n=1:hex_size
    model.geom('geom1').feature().create(strcat('h',num2str(n)),'Hexahedron');
    model.geom('geom1').feature(strcat('h',num2str(n))).set('p',[
physical_hex_info{n}(3,1),physical_hex_info{n}(7,1),physical_hex_info{n}(8,1),physical_hex_info{n}(4,1),physical_hex_info{n}(1,1),physical_hex_info{n}(5,1),physical_hex_info{n}(6,1),physical_hex_info{n}(2,1);
physical_hex_info{n}(3,3),physical_hex_info{n}(7,3),physical_hex_info{n}(8,3),physical_hex_info{n}(4,3),physical_hex_info{n}(1,3),physical_hex_info{n}(5,3),physical_hex_info{n}(6,3),physical_hex_info{n}(2,3);
physical_hex_info{n}(3,2),physical_hex_info{n}(7,2),physical_hex_info{n}(8,2),physical_hex_info{n}(4,2),physical_hex_info{n}(1,2),physical_hex_info{n}(5,2),physical_hex_info{n}(6,2),physical_hex_info{n}(2,2)]);
end
%Unioning all of the subdomains to reduce the "object creation overhead"
model.geom('geom1').feature().create('h1 h2','Union');

model.geom('geom1').run();
mphgeom(model,'geom1')




%Mesh is Currently Disabled due to Java Heap space Shortage

model.geom('geom1').run();
% mesh = model.mesh.create('mesh','geom1');
% size = mesh.feature('size');
% size.set('hmax','mh');
% size.set('hmin','mh-mh/10');
% size.set('hcurve','0.2');
% ftet = mesh.feature.create('ftet','FreeTet');
% mesh.run;

% figure
% mphmesh(model)
 
ModelUtil.showProgress(true);

t_comsol = toc;
disp(sprintf('t_comsol: %d sec.',t_comsol))

% num_all_elements = model.mesh('mesh').stat().getNumElem('all');
% disp(sprintf('Number of all elements: %d', num_all_elements));

% x=1:length(a_physical_form);
% figure
% plot(x,a1_virtual,'r',x,a2_virtual,'b',x,a3_virtual,'g',x,a4_virtual,'k')
% title(sprintf('a_n Edge Lengths (%dx%dx%d Volume Mesh)',i-1,j-1,k-1),'FontSize',12,'Fontweight','bold')
% xlabel('Edge index','FontSize',11,'Fontweight','bold')
% ylabel('Edge Lengths (mm)','FontSize',11,'Fontweight','bold')
% xlim([1 length(a_physical_form)])
% legend('a_1','a_2','a_3','a_4')
% 
% x=1:length(b_physical_form);
% figure
% plot(x,b1_virtual,'r',x,b2_virtual,'b',x,b3_virtual,'g',x,b4_virtual,'k')
% title(sprintf('b_n Edge Lengths (%dx%dx%d Volume Mesh)',i-1,j-1,k-1),'FontSize',12,'Fontweight','bold')
% xlabel('Edge index','FontSize',11,'Fontweight','bold')
% ylabel('Edge Lengths (mm)','FontSize',11,'Fontweight','bold')
% xlim([1 length(a_physical_form)])
% legend('b_1','b_2','b_3','b_4')
% 
% x=1:length(c_physical_form);
% figure
% plot(x,c1_virtual,'r',x,c2_virtual,'b',x,c3_virtual,'g',x,c4_virtual,'k')
% title(sprintf('c_n Edge Lengths (%dx%dx%d Volume Mesh)',i-1,j-1,k-1),'FontSize',12,'Fontweight','bold')
% xlabel('Edge index','FontSize',11,'Fontweight','bold')
% ylabel('Edge Lengths (mm)','FontSize',11,'Fontweight','bold')
% xlim([1 length(a_physical_form)])
% legend('c_1','c_2','c_3','c_4')
% 
% x=1:length(cuboid_volume1_virtual);
% figure
% plot(x,cuboid_volume1_virtual,'r',x,cuboid_volume2_virtual,'b',x,cuboid_volume3_virtual,'g',x,hex_volume_virtual,'k')
% title(sprintf('Block Volumes (%dx%dx%d Volume Mesh) Virtual',i-1,j-1,k-1),'FontSize',12,'Fontweight','bold')
% xlabel('Block index','FontSize',11,'Fontweight','bold')
% ylabel('Volume (mm^3)','FontSize',11,'Fontweight','bold')
% xlim([1 length(a_physical_form)])
% legend('Cuboid-1','Cuboid-2','Cuboid-3','Hexahedral')
% 
% 
% x=1:length(a_physical_form);
% figure
% plot(x,a1_physical,'r',x,a2_physical,'b',x,a3_physical,'g',x,a4_physical,'k')
% title(sprintf('a_n Edge Lengths (%dx%dx%d Volume Mesh)',i-1,j-1,k-1),'FontSize',12,'Fontweight','bold')
% xlabel('Edge index','FontSize',11,'Fontweight','bold')
% ylabel('Edge Lengths (mm)','FontSize',11,'Fontweight','bold')
% xlim([1 length(a_physical_form)])
% legend('a_1','a_2','a_3','a_4')
% 
% x=1:length(b_physical_form);
% figure
% plot(x,b1_physical,'r',x,b2_physical,'b',x,b3_physical,'g',x,b4_physical,'k')
% title(sprintf('b_n Edge Lengths (%dx%dx%d Volume Mesh)',i-1,j-1,k-1),'FontSize',12,'Fontweight','bold')
% xlabel('Edge index','FontSize',11,'Fontweight','bold')
% ylabel('Edge Lengths (mm)','FontSize',11,'Fontweight','bold')
% xlim([1 length(a_physical_form)])
% legend('b_1','b_2','b_3','b_4')
% 
% x=1:length(c_physical_form);
% figure
% plot(x,c1_physical,'r',x,c2_physical,'b',x,c3_physical,'g',x,c4_physical,'k')
% title(sprintf('c_n Edge Lengths (%dx%dx%d Volume Mesh)',i-1,j-1,k-1),'FontSize',12,'Fontweight','bold')
% xlabel('Edge index','FontSize',11,'Fontweight','bold')
% ylabel('Edge Lengths (mm)','FontSize',11,'Fontweight','bold')
% xlim([1 length(a_physical_form)])
% legend('c_1','c_2','c_3','c_4')
% 
% x=1:length(cuboid_volume1_physical);
% figure
% plot(x,cuboid_volume1_physical,'r',x,cuboid_volume2_physical,'b',x,cuboid_volume3_physical,'g',x,hex_volume_physical,'k')
% title(sprintf('Block Volumes (%dx%dx%d Volume Mesh)',i-1,j-1,k-1),'FontSize',12,'Fontweight','bold')
% xlabel('Block index','FontSize',11,'Fontweight','bold')
% ylabel('Volume (mm^3)','FontSize',11,'Fontweight','bold')
% xlim([1 length(a_physical_form)])
% legend('Cuboid-1','Cuboid-2','Cuboid-3','Hexahedral')



