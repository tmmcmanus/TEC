function [gamma_111,gamma_112,gamma_221,gamma_222,gamma_121,gamma_122] = Chris_first_kind(pdxh11,pdxh12,pdyh12,pdyh22,pdyh11,pdxh22);

gamma_111 = 0.5*(pdxh11);
gamma_112 = 0.5*(2*(pdxh12) - pdyh11);
gamma_221 = 0.5*(2*(pdyh12) - pdxh22);
gamma_222 = 0.5*(pdyh22);
gamma_121 = 0.5*(pdyh11);
gamma_122 = 0.5*(pdxh22);
end