% EXPERIMENT 1

dir1 = ['data/experiments/', 'fulldimNrscan', '/', 'parallel', '/'];
[errs1,nums1,pids1,labels1] = gather_plotdata(dir1, 32);

tr_range = 1;
dir2 = ['data/experiments/', 'fulldimNr10000-12000', '/', 'parallel', '/'];
[errs2,nums2,pids2,labels2] = gather_plotdata(dir2, 8, tr_range);

dir3 = ['data/experiments/', 'fulldimNr12000', '/', 'parallel', '/'];
[errs3,nums3,pids3,labels3] = gather_plotdata(dir3, 8);

figure(1)
my_boxplot([nums1, nums2, nums3]);
xticklabels({'2000','4000','6000','8000', ...
             num2str(labels2.hyp_range(1)), ...
             num2str(labels3.hyp_range(1))})

xlabel(labels2.xlab)
ylabel(labels2.ylab)

description = sprintf('Full dimensional \n hybrid setup');
text(1.2*min(xlim),0.9*max(ylim),description,'fontsize',14)

figure(2)
E = [errs1, errs2, errs3];
for i = 1:1:numel(E)
    plot(E{i}(1:end-1))
    hold on
end
plot(xlim, [4,4])
hold off
xlabel('Days')
ylabel('Error')
ylim([0,4.5])
title('Error evolution')

dat62 = load([dir2, '/results_5.mat']);
figure(3)
nx  = 64; ny = 64; nun = 2;
dim = nx*ny*nun;
Ldim    = 1e6;
Udim    = 3.171e-2;
tdim    = Ldim / Udim; % in seconds
scaling = 3600*24/tdim;
subplot(1,2,1)
plotQG(nx,ny,1,scaling*dat62.predictions{62},false)
title('predicted vorticity (days^{-1})')
subplot(1,2,2)
plotQG(nx,ny,1,scaling*dat62.truths{62},false)
title('true vorticity (days^{-1})')

figure(4)
qg = QG(nx, ny, 1 );  % coarse QG with periodic bdc
qg.set_par(18, 0  );  % stirring type: 0 = cos(5x), 1 = sin(16x)
qg.set_par(11, 0.5);  % stirring amplitude
qg.set_par(5,  100);  % Reynolds number for coarse model
cols = lines(10);
f1 = plotQGspectrum(qg, scaling*dat62.truths{62}, 5, cols(1,:));
hold on;
f2 = plotQGspectrum(qg, scaling*dat62.predictions{62}, 5, cols(2,:));
hold off;