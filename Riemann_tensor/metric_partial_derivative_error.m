function[X_f2,Y_f2,pdxh11_true,pdxh12_true,pdxh22_true,pdyh11_true,pdyh12_true,pdyh22_true,pdxh11_error,pdxh12_error,pdxh22_error,pdyh11_error,pdyh12_error,pdyh22_error] = metric_partial_derivative_error(X_f1,Y_f1,i,pdxh11,pdxh12,pdxh22,pdyh11,pdyh12,pdyh22);

X_f2 = X_f1;
X_f2(1:i-2)=[];
X_f2(1:i-2:end)=[];
X_f2(i-3:i-3:end)=[];
X_f2(end-(i-5):end)=[];

Y_f2 = Y_f1;
Y_f2(1:i-2)=[];
Y_f2(1:i-2:end)=[];
Y_f2(i-3:i-3:end)=[];
Y_f2(end-(i-5):end)=[];

pdxh11_true = zeros(length(Y_f2),1);
pdxh12_true = Y_f2;
pdxh22_true = 2.*(X_f2);
pdyh11_true = 2.*(Y_f2);
pdyh12_true = X_f2;
pdyh22_true = zeros(length(X_f2),1);


pdxh11_numeric = pdxh11;
pdxh12_numeric = pdxh12;
pdxh22_numeric = pdxh22;
pdyh11_numeric = pdyh11;
pdyh12_numeric = pdyh12;
pdyh22_numeric = pdyh22;

pdxh11_error = abs(pdxh11_true - pdxh11');
pdxh12_error = abs(pdxh12_true - pdxh12');
pdxh22_error = abs(pdxh22_true - pdxh22');
pdyh11_error = abs(pdyh11_true - pdyh11');
pdyh12_error = abs(pdyh12_true - pdyh12');
pdyh22_error = abs(pdyh22_true - pdyh22');
end


