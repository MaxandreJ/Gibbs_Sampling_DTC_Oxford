K= 12;
n_iterations= 500;
burn_in=100;
a = 0.1;
mu_start = 0.5;
mu_unknown = 1;
beta= [1,1];
[ Z, S, mu, max_lr, min_ent, min_ent_M, min_ent_s, max_lr_M,max_lr_s, posterior_mean_M, information,background ]  = find_motifs('data/data3.fasta',K, n_iterations,burn_in, a, mu_start, mu_unknown, beta)
plot(mu)
xlabel('iteration')
ylabel('\mu')
spy(Z)

%%task4
[ Z, S, mu, max_lr, min_ent, min_ent_M, min_ent_s, max_lr_M,max_lr_s, posterior_mean_M, information,background ]  = find_motifs_phylo('data/data4.fasta',11, 500,100, 0.2, 1, 0, [2,2])
