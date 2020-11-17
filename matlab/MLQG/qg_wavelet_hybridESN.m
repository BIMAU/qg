compute = true;
enslave = false;
hybrid  = false;
movie   = false;

Ldim    = 1e6;
Udim    = 3.171e-2;
tdim    = Ldim / Udim; % in seconds
scaling = 3600*24/tdim;

if compute
    fprintf('load training data...\n'); tic;
    % fname_base = 'N256-N128_ff2_Re4.0e+04-Re4.0e+02_Tstart151_Tend167';
    fname_base = 'N128-N64_ff2_Re1.0e+04-Re1.0e+02_Tstart159_Tend187';
    % trdata = load(['data/training/', fname_base, '.mat']);
    fprintf('load training data... done (%fs)\n', toc);

    nxc  = trdata.nxc;
    nyc  = trdata.nyc;
    nun  = 2;
    dim  = nxc*nyc*nun;
    RX   = trdata.RX;
    ERX  = trdata.ERX;
    ampl = trdata.ampl;
    stir = trdata.stir;
    Re_c = trdata.Re_c;

    %%
    % create coarse QG setup for this problem
    qgc = QG(nxc, nyc, 1);  % coarse QG with periodic bdc
    qgc.set_par(18, stir);  % stirring type: 0 = cos(5x), 1 = sin(16x)
    qgc.set_par(11, ampl);  % stirring amplitude
    qgc.set_par(5,  Re_c);  % Reynolds number for coarse model

    % dimension reduction Na
    Na = dim / 2;

    % create wavelet basis, set block size
    bs = 32;
    H  = create_wavelet_basis(nxc, nyc, nun, bs, true);
    Ha = H(1:Na,:)';
    Hd = H(Na+1:dim,:)';
    Ph = Ha*Ha';

    % Alternative: POD basis
    % [U,~,~] = svd(RX);
    % Ha = U(1:Na,:)';
    % Hd = U(Na+1:dim,:)';

    % create input/output data for hybrid ESN
    if hybrid
        U = [Ha' * RX(:,1:end-1); Ha' * ERX(:,1:end-1)];
    else
        U = [Ha' * RX(:,1:end-1)];
    end

    Y = [Ha' * RX(:,2:end )];

    % separate training and testing data
    cutoff = round(trdata.tpars.Nt * 3 / 4)
    fprintf(' training time: %d\n', cutoff)
    fprintf(' test time: %d\n', trdata.tpars.Nt - cutoff)
    trainU = U(:,1:cutoff)';
    trainY = Y(:,1:cutoff)';
    testU  = U(:,cutoff+1:end)';
    fullY  = RX(:,2:end);
    fullY  = fullY(:,cutoff+1:end);
    testY  = fullY';

    % ESN parameters
    pars                    = {};
    pars.scalingType        = 'standardize';
    pars.Nr                 = 3000;
    pars.rhoMax             = 0.3;
    pars.alpha              = 1.0;
    pars.Wconstruction      = 'avgDegree';
    pars.avgDegree          = 10;
    pars.lambda             = 1e-1;
    pars.bias               = 0.0;
    pars.squaredStates      = 'even';
    pars.reservoirStateInit = 'random';
    pars.inputMatrixType    = 'balancedSparse';
    pars.inAmplitude        = 1.0;

    if hybrid
        pars.feedThrough        = true;
        pars.ftRange            = Na+1:2*Na;
    end

    esn = ESN(pars.Nr, size(U,1), size(Y,1));
    esn.setPars(pars);
    esn.initialize;
    esn.train(trainU, trainY);

    Npred = size(testU, 1);
    Npred = 80;
    state = esn.X(end,:);
    predY = zeros(Npred, dim);
    predS = zeros(Npred, dim);
    [~, ~, nSpec, ~] = computeQGspectrum(qgc, zeros(dim,1));
    testSpec  = zeros(Npred, nSpec);
    predSpecY = zeros(Npred, nSpec);
    predSpecS = zeros(Npred, nSpec);
    err_hyb = zeros(Npred,1);
    err_std = zeros(Npred,1);
    yk  = RX(:,cutoff+1);
    yks = RX(:,cutoff+1);

    for i = 1:Npred
        fprintf('step %d/%d\n', i, Npred);
        % coarse model prediction:
        Eyk   = qgc.step(yk, trdata.tpars.dt);

        % append prediction to regular input data
        if hybrid
            u_in  = [Ha'*yk(:); Ha'*Eyk(:)]';
        else
            u_in  = [Ha'*yk(:)]';
        end
        u_in  = esn.scaleInput(u_in);
        state = esn.update(state, u_in)';
        u_out = esn.apply(state, u_in);
        u_out = esn.unscaleOutput(u_out);

        if enslave
            cd = Hd'*Eyk;
            for k = 1:1
                yk = Ha*u_out(:) + Hd*cd;
                Gd = Hd'*qgc.rhs(yk);
                Jd = Hd'*qgc.jacobian(yk, 0.0)*Hd;
                up = Jd \ Gd;
                cd = cd + up;
                fprintf('Projected iteration %e\n',norm(up,2));
                if norm(up,2) < 1e-2
                    break;
                end
            end

            yk = Ha*u_out(:) + Hd*cd;
        else
            % kappa = 1-min(mod(i,2),1)
            kappa = 1;
            yk = kappa*Ha*u_out(:) + Hd*(Hd'*Eyk) + (1-kappa)*Ha*(Ha'*Eyk);
        end

        predY(i,:) = yk;

        % standalone step ---
        yks = qgc.step(yks, trdata.tpars.dt);
        predS(i,:) = yks;
        
        % compute spectra
        testSpec(i,:)  = computeQGspectrum(qgc, scaling*testY(i,:));
        predSpecY(i,:) = computeQGspectrum(qgc, scaling*predY(i,:));
        predSpecS(i,:) = computeQGspectrum(qgc, scaling*predS(i,:));

        % decide where to test in the spectrum
        sd = find(testSpec(i,:) > 1e-10);

        diff_hyb   = (log(testSpec(i,sd)) - log(predSpecY(i,sd)));
        diff_std   = (log(testSpec(i,sd)) - log(predSpecS(i,sd)));
        err_hyb(i) = norm(diff_hyb, 2);
        err_std(i) = norm(diff_std, 2);
    end

    % plotting -------------
    subplot(2,3,1)
    plotQG(nxc, nyc, 1, scaling*fullY(:,i), false)

    subplot(2,3,3)

    hold on
    plot(err_hyb(1:i)); hold on
    plot(err_std(1:i)); 
    plot(xlim, [5, 5]);  % ARBITRARY #FIXME
    hold off

    subplot(2,3,6)
    cols = lines(10);
    f1 = plotQGspectrum(qgc, nxc, nyc, scaling*fullY(:,i), 5, cols(1,:));
    hold on;
    f2 = plotQGspectrum(qgc, nxc, nyc, scaling*predS(i,:), 5, cols(2,:));
    hold on;
    f3 = plotQGspectrum(qgc, nxc, nyc, scaling*predY(i,:), 5, cols(3,:));
    hold off;
    if hybrid
        legend([f1,f2,f3],'data', 'coarse', 'hybrid');
    else
        legend([f1,f2,f3],'data', 'coarse', 'ESN');
    end
    ylim([10^(-10), 10^5]);

    subplot(2,3,4)
    plotQG(nxc,nyc,1,scaling*yks,false)

    subplot(2,3,5)
    plotQG(nxc,nyc,1,scaling*yk,false)
end

if movie
    if hybrid
        fname = sprintf('data/avi/N64qg_wrom%d_W%d_hybrid.avi',...
                        dim/Na,pars.Nr);
    else
        fname = sprintf('data/avi/N64qg_wrom%d_W%d_standl.avi',...
                        dim/Na,pars.Nr);
    end

    fprintf('creating movie: %s\n',fname);
    writerObj = VideoWriter(fname, 'Motion JPEG AVI');
    writerObj.FrameRate = 16;
    writerObj.Quality = 90;
    open(writerObj);
    set(0,'DefaultFigureWindowStyle','normal')
    fhandle = figure('units','pixels','position',[100,100,1500,800]);
    set(gca,'position',[0.05 0.1 .92 0.85],'units','normalized');
    set(gca,'color','w','fontsize',15);
else
    return
end


for i = 1:2:Npred
    subplot(2,3,1)
    plotQG(nxc,nyc,1,scaling*fullY(:,i),false)
    caxis([-0.2,0.2])
    title(['restricted data, day: ', ...
           num2str(i)])

    subplot(2,3,4)
    [u,v]  = qgc.compute_uv(scaling*fullY(:,i));
    u      = reshape(u,nxc,nyc);
    v      = reshape(v,nxc,nyc);
    Efield = (u.^2+v.^2)';
    imagesc(Efield)
    set(gca,'ydir','normal')
    titleString = sprintf('u^2 + v^2, sum = %1.2e', sum(sum(Efield)));
    title(titleString);
    colorbar
    caxis([0,12e-5]);

    subplot(2,3,5)
    [u,v]  = qgc.compute_uv(scaling*predY(i,:));
    u      = reshape(u,nxc,nyc);
    v      = reshape(v,nxc,nyc);
    Efield = (u.^2+v.^2)';
    imagesc(Efield)
    set(gca,'ydir','normal')
    titleString = sprintf('u^2 + v^2, sum = %1.2e', sum(sum(Efield)));
    title(titleString);
    colorbar
    caxis([0,12e-5]);

    subplot(2,3,3)
    cols = lines(10);
    f1 = plotQGspectrum(qgc, nxc, nyc, scaling*fullY(:,i), 5, cols(1,:));
    hold on;
    f2 = plotQGspectrum(qgc, nxc, nyc, scaling*predS(i,:), 5, cols(2,:));
    hold on;
    f3 = plotQGspectrum(qgc, nxc, nyc, scaling*predY(i,:), 5, cols(3,:));
    hold off;
    if hybrid
        legend([f1,f2,f3],'data', 'coarse', 'hybrid');
    else
        legend([f1,f2,f3],'data', 'coarse', 'ESN');
    end
    ylim([10^(-10), 10^5]);

    subplot(2,3,2)
    plotQG(nxc,nyc,1,scaling*predY(i,:),false)
    caxis([-0.2,0.2])
    if hybrid
        title(['hybrid ESN + wavelet reduction: ', num2str(dim/Na)])
    else
        title(['ESN + wavelet reduction: ', num2str(dim/Na)])
    end

    subplot(2,3,6)
    plotQG(nxc,nyc,1,scaling*predS(i,:),false)
    caxis([-0.2,0.2])
    title(['standalone coarse model'])

    if movie
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    else
        drawnow
    end
end

if movie
    close(writerObj);
end