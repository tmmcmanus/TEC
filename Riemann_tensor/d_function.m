function [d,d_vector] = d_function(X,Y,Z,i,j)
    for n=1:(i.*(j-1))
    d(n)=sqrt((X(n) - X(n+i)).^2 +(Y(n)-Y(n+i)).^2+(Z(n)-Z(n+i)).^2);
    d_vector(n,:)=[(X(n)-X(n+i)),(Y(n)-Y(n+i)),(Z(n)-Z(n+i))];
    end
end