function[gamma_111_error,gamma_112_error,gamma_221_error,gamma_222_error,gamma_121_error,gamma_122_error,gamma_111_true,gamma_112_true,gamma_221_true,gamma_222_true,gamma_121_true,gamma_122_true] = Chris_first_kind_error(pdxh12_true, pdxh11_true, pdyh12_true, pdxh22_true,pdyh22_true,pdyh11_true,gamma_111,gamma_112,gamma_221,gamma_222,gamma_121,gamma_122);
gamma_111_true = 0.5*(pdxh11_true);
gamma_112_true = 0.5*(2*(pdxh12_true) - pdyh11_true);
gamma_221_true = 0.5*(2*(pdyh12_true) - pdxh22_true);
gamma_222_true = 0.5*(pdyh22_true);
gamma_121_true = 0.5*(pdyh11_true);
gamma_122_true = 0.5*(pdxh22_true);

gamma_111_error = abs(gamma_111_true - gamma_111');
gamma_112_error = abs(gamma_112_true - gamma_112');
gamma_221_error = abs(gamma_221_true - gamma_221');
gamma_222_error = abs(gamma_222_true - gamma_222');
gamma_121_error = abs(gamma_121_true - gamma_121');
gamma_122_error = abs(gamma_122_true - gamma_122');
end