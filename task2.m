max_lrs= [];
min_ents=[];
min_ent_Ms = {};
max_lr_Ms = {};
posterior_mean_Ms = {};
% for K=5:15
for a=[0.01,0.1]
  [ Z, S, mu, max_lr, min_ent, ...
             min_ent_M, min_ent_s, ...
             max_lr_M,max_lr_s, ...
             posterior_mean_M, information, bacgkround ]  = find_motifs('data2.fasta',12, ...
                                                        200,100, ...
                                                        a, 1, 0);

  max_lrs = [max_lrs max_lr*(1/4)^(K+K-5)];
  min_ents = [min_ents min_ent];
  min_ent_Ms =[min_ent_Ms; {min_ent_M}];
  max_lr_Ms  = [max_lr_Ms; {max_lr_M}];
  posterior_mean_Ms  = [posterior_mean_Ms; {posterior_mean_M}];
end
