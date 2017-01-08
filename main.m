%%% Script to run the Gibbs Sampler (instead of copy/pasting the function
%%% all the time...)

% Set parameters
sequence_file= 'data/data1.fasta';
K= 11;
n_iterations= 200;
burn_in= 100;
a = 1;
mu_start = 1;
mu_unknown = 1;
beta= [1,1];

% Find the motifs by running the Gibbs sampler
[ Z, S, mu, max_lr, min_ent, ...
           min_ent_M, min_ent_s, ...
           max_lr_M,max_lr_s, ...
           posterior_mean_M, information,background ]  = find_motifs(sequence_file,K, ...
                                                      n_iterations,burn_in, ...
                                                      a, mu_start, mu_unknown, beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES ON ARGUMENTS
%
%% sequence_file: a FASTA-formatted file containing the input sequences
%% K:             the length of the motif
%% n_iterations:  number of iterations for which Gibbs sampler
%                 should be run
%% burn_in:       number of iterations to allow for the burn-in
%                 phase, while the MCMC is converging.
%% a:             a constant multiplier for the uniform prior on the
%                 motif
%% mu_start:      the starting value of mu
%% mu_unknown:    0 if mu is fixed to mu_start, 1 if mu is unknown and
%                 is to be estimated by MCMC.
%% beta:          a length-two vector, containing prior parameters
%                 for mu (used only if mu_unknown == 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Display the results
% Display the M for which the likelihood ratio is maximal
seqlogo_fig1(max_lr_M)
% Display the average information per site as a function of iterations
plot(information)
ylabel('Average information per site')
xlabel('Iteration')
