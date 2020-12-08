[errs1,nums1,pids1,labels1] = gather_plotdata('fulldimNrscan', ...
                                              'parallel', 32);

tr_range = 1;
[errs2,nums2,pids2,labels2] = gather_plotdata('fulldimNr10000-12000', ...
                                              'parallel', 8, tr_range);


my_boxplot([nums1, nums2]);
xticklabels({'2000','4000','6000','8000', ...
             num2str(labels2.hyp_range(1))})
xlabel(labels2.xlab)
ylabel(labels2.ylab)