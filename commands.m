sequence_file= 'data/data2.fasta';
K= 11;
n_iterations= 10000;
burn_in= 2000;
% a = 0.1;
mu_start = 1;
mu_unknown = 0;
beta= [1,1];

for a=[0.1]

%%task4
% [ Z, S, mu, max_lr, min_ent, min_ent_M, min_ent_s, max_lr_M,max_lr_s, posterior_mean_M, information,background ]  = find_motifs_phylo(sequence_file,K,n_iterations,burn_in, a, mu_start, mu_unknown, beta)

%%task5
[ Z, S, mu, max_lr, min_ent, min_ent_M, min_ent_s, max_lr_M,max_lr_s, posterior_mean_M, information,background ]  = find_motifs_K(sequence_file,K,n_iterations,burn_in, a, mu_start, mu_unknown, beta)


figure();
plot(information)
title(['Parameter values: a =', num2str(a),', K =', num2str(K),...
    ', burn\_in =', num2str(burn_in), ', n\_iterations =', num2str(n_iterations), ...
    ', mu\_start =', num2str(mu_start), ', mu\_unknown =', num2str(mu_unknown), ...
    ', sequence\_file =', sequence_file, ', beta =', num2str(beta)]);
ylabel('Average information per site')
xlabel('Iteration')

[figure_handle,~] = seqlogo_fig1(max_lr_M);
figure_handle.Name = ['a =', num2str(a),', K =', num2str(K),...
    ', burn_in =', num2str(burn_in), ', n_iterations =', num2str(n_iterations), ...
    ', mu_start =', num2str(mu_start), ', mu_unknown =', num2str(mu_unknown), ...
    ', sequence_file =', sequence_file, ', beta =', num2str(beta)];

end
