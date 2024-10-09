% using the function: cosmo_balance_partitions from cosmoMvpa


function acc = searchlight_func(ds,measure_args)
acc = struct();
ds=cosmo_remove_useless_data(ds);

%% set the targets and chunks
ds.sa.targets = [ones(measure_args.s1_H,1);ones(measure_args.s1_L,1)+1;ones(measure_args.s2_H,1);ones(measure_args.s2_L,1)+1];     %train里面20来个1，2； test组里也是20来个1，2 也可能30个全了
ds.sa.chunks = [ones(measure_args.s1_H + measure_args.s1_L,1);ones(measure_args.s2_H+measure_args.s2_L,1)+1];  % 应该学习阶段的可用试次都是chunk1，自创阶段的都是chunk2

x=cosmo_nfold_partitioner(ds);
aa =cosmo_balance_partitions(x,ds,'nmin',1);

counter=1;

for m=1:length(aa.train_indices)
    
   if aa.train_indices{1,m}(1,1)==1
    train=cosmo_slice(ds,aa.train_indices{1,m}',1);
    test=cosmo_slice(ds,aa.test_indices{1,m}',1);
%%  
    pred = cosmo_classify_svm(train.samples, train.sa.targets, test.samples);
    acc0 = mean(test.sa.targets == pred);
    myacc(counter,1)=acc0;
    counter=counter+1;
   end
    
end

acc.samples = mean(myacc)-0.5;
end

