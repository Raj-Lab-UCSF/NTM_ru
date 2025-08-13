function [] = NTM_sim(matdir, seed_in, study_in)

% matdir = '/Users/nbarron/Desktop/NTM_ru/Transport_Network_Neumann/MatFiles';

rng(seed_in);

r_alpha = rand();
r_gamma = rand();
r_delta = rand();
r_epsilon = rand();
r_uprel = rand();

alpha_in = r_alpha; % 0 - 1
gamma_in = ((0.008 - 0.001) * r_gamma) + 0.001; % 0.001 - 0.008
delta_in = ((100 - 10) * r_delta) + 10; % 10 - 100
epsilon_in = ((100 - 10) * r_epsilon) + 10; % 10 - 100
uprel_in = ((778 - 50) * r_uprel) + 50; % 50 - 778

save_file = "/Users/nbarron/Desktop/simulation_" + num2str(seed_in) + ".mat";

call_NTM2(matdir, save_file, study_in, alpha_in, gamma_in, delta_in, epsilon_in, uprel_in);