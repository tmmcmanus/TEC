function [a,a_vector] = a_function(X,Y,Z,i,j)
for n = 1:(i.*j)-1
    a(n)= sqrt((X(n)-X(n+1)).^2 + (Y(n)-Y(n+1)).^2 + (Z(n) - Z(n+1)).^2);
    a_vector(n,:)=[(X(n)-X(n+1)),(Y(n)-Y(n+1)),(Z(n) - Z(n+1))];
end
a(i:i:end)=[]; %removing false lengths that are a result of populating scheme 
a_vector(i:i:end,:)=[];
end