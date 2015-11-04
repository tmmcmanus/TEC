function[gamma_1_11_error,gamma_1_21_error,gamma_1_22_error,gamma_2_11_error,gamma_2_12_error,gamma_2_22_error] = Chris_second_kind_error(gamma_1_11,gamma_1_21,gamma_1_22,gamma_2_11,gamma_2_12,gamma_2_22,gamma_1_11_true,gamma_1_21_true,gamma_1_22_true,gamma_2_11_true,gamma_2_12_true,gamma_2_22_true);

gamma_1_11_error = abs(gamma_1_11_true - gamma_1_11');
gamma_1_21_error = abs(gamma_1_21_true - gamma_1_21');
gamma_1_22_error = abs(gamma_1_22_true - gamma_1_22');
gamma_2_11_error = abs(gamma_2_11_true - gamma_2_11');
gamma_2_12_error = abs(gamma_2_12_true - gamma_2_12');
gamma_2_22_error = abs(gamma_2_22_true - gamma_2_22');
end