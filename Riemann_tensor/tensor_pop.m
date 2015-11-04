function[covariant, contravariant] = tensor_pop(h11_central,h12_central,h22_central,i,j)
covariant = cell(i-2);
contravariant = cell(i-2);
h11_central_matrix = reshape(h11_central,i-2,j-2);
h12_central_matrix = reshape(h12_central,i-2,j-2);
h22_central_matrix = reshape(h22_central,i-2,j-2);

for m = 1:i-2
    for n = 1:j-2
      covariant{m,n}=[h11_central_matrix(m,n),h12_central_matrix(m,n);h12_central_matrix(m,n),h22_central_matrix(m,n)];
      contravariant{m,n}=inv([h11_central_matrix(m,n),h12_central_matrix(m,n);h12_central_matrix(m,n),h22_central_matrix(m,n)]);
    end
end      
end