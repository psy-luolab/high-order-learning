function my_rd=searchlight_Func(my_ds,args)
my_rd=struct();
total_explained_hsc=0;
my_dsm=cosmo_pdist(my_ds.samples, 'correlation');
RDM=cosmo_squareform(my_dsm);

[~,~,explained]=pcacov(RDM);

for n = 1:length(explained)
    total_explained_hsc=total_explained_hsc+explained(n);
    if total_explained_hsc > 70
        break
    end
end

RD_var_hsc = n;

my_rd.samples=RD_var_hsc;
end
