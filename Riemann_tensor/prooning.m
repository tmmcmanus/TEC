function [a,a_vector,c,c_vector,b,b_vector,d,d_vector] = prooning(a,a_vector,c,c_vector,b,b_vector,d,d_vector,i,j)
a(((i-1).*(j-1))+1:((i-1).*j))=[];
a_vector(((i-1).*(j-1))+1:((i-1).*j),:)=[];
c(1:(i-1))=[];
c_vector(1:(i-1),:)=[];
b(1:i:end)=[];
b_vector(1:i:end,:)=[];
d(i:i:end)=[];
d_vector(i:i:end,:)=[];
end
