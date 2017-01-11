%%% Script to run the Gibbs Sampler (instead of copy/pasting the function
%%% all the time...)

% Set parameters

sequence_file= 'data/data3.fasta';
%K= 20 ;
n_iterations= 2000;
burn_in= 1000;
a = 0.1;
mu_start = 1;
mu_unknown = 1;
beta= [2,4];


% Find the motifs by running the Gibbs sampler

max_lrs=[];
min_ents=[];
min_ent_Ms = {};
max_lr_Ms = {};
posterior_mean_Ms = {};
informations = {};

for mu_start = 0.8%[0.3, 0.5, 0.8, 1]
for K=[18]

% for K=5:15

for a=[0.05, 0.1, 0.2, 1, 2, 15]


  [ Z, S, mu, max_lr, min_ent, ...
             min_ent_M, min_ent_s, ...
             max_lr_M,max_lr_s, ...
             posterior_mean_M, information,background ]  = find_motifs(sequence_file,K, ...
                                                        n_iterations,burn_in, ...
                                                        a, mu_start, mu_unknown, beta)

  max_lrs = [max_lrs max_lr*(1/4)^(K+K-5)];
  min_ents = [min_ents min_ent];
  min_ent_Ms =[min_ent_Ms; {min_ent_M}];
  max_lr_Ms  = [max_lr_Ms; {max_lr_M}];
  posterior_mean_Ms  = [posterior_mean_Ms; {posterior_mean_M}];
  informations  = [informations; {information}];


  % figure();,
  [figure_handle,~] = seqlogo_fig1(posterior_mean_M);
  figure_handle.Name = ['a =', num2str(a),', K =', num2str(K),...
      ', burn_in =', num2str(burn_in), ', n_iterations =', num2str(n_iterations), ...
      ', mu_start =', num2str(mu_start), ', mu_unknown =', num2str(mu_unknown), ...
      ', sequence_file =', sequence_file, ', beta =', num2str(beta)];


  % Display the average information per site as a function of iterations
  figure();
  plot(information)
  title(['Parameter values: a =', num2str(a),', K =', num2str(K),...
      ', burn\_in =', num2str(burn_in), ', n\_iterations =', num2str(n_iterations), ...
      ', mu\_start =', num2str(mu_start), ', mu\_unknown =', num2str(mu_unknown), ...
      ', sequence\_file =', sequence_file, ', beta =', num2str(beta)]);
  ylabel('Average information per site')
  xlabel('Iteration')
  
  
    figure();
bar(mean(Z(:,burn_in:end),2))
  title(['Parameter values: a =', num2str(a),', K =', num2str(K),...
      ', burn\_in =', num2str(burn_in), ', n\_iterations =', num2str(n_iterations), ...
      ', mu\_start =', num2str(mu_start), ', mu\_unknown =', num2str(mu_unknown), ...
      ', sequence\_file =', sequence_file, ', beta =', num2str(beta)]);
  ylabel('Probability that the sequence is selected')
  xlabel('Sequence number')
  
  % plot(mu)
  % xlabel('iteration')
  % ylabel('\mu')
  % spy(Z)
end
end
end


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
% Display2 the M for which the likelihood ratio is maximal

% seqlogo_fig1(max_lr_M)
% % Display the average information per site as a function of iterations
% plot(information)
% ylabel('Average information per site')
% xlabel('Iteration')

%% Task 1. done
%%% Sequence logo graph for mean and max likelihood.
%%% convergence
%%% a sweep
%%%% plot convergences, and lr, ents.

%% Task 2. See data2results.
%%%
%% Task 3. Learning mu.

%% Task 4. Phylogenetic information.

%% Task 5.
