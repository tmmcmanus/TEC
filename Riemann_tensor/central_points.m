function [X_central,Y_central,Z_central]=central_points(X,Y,Z,i,j)
X_central=X; Y_central=Y; Z_central=Z;
X_central((i.*(j-1))+1:(i.*j))=[];
X_central(i:i:end)=[];
X_central(1:i-1:end)=[];
X_central(1:i-2)=[];
Y_central((i.*(j-1))+1:(i.*j))=[];
Y_central(i:i:end)=[];
Y_central(1:i-1:end)=[];
Y_central(1:i-2)=[];
Z_central((i.*(j-1))+1:(i.*j))=[];
Z_central(i:i:end)=[];
Z_central(1:i-1:end)=[];
Z_central(1:i-2)=[];
end
