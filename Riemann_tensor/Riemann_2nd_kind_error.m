function[R_2_112_error, R_1_112_error,R_2_221_error,R_1_212_error] = Riemann_2nd_kind_error(R_2_112_true,R_1_112_true,R_2_221_true,R_1_212_true,R_2_112,R_1_112,R_2_221,R_1_212);
R_2_112_error = abs(R_2_112_true - R_2_112');
R_1_112_error = abs(R_1_112_true - R_1_112');
R_2_221_error = abs(R_2_221_true - R_2_221');
R_1_212_error = abs(R_1_212_true - R_1_212');
end
