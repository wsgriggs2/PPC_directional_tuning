function [imageOut,template1] = normcorre_doppler(imageIn, verbose,template)
gcp;

if ~exist('verbose','var')
    verbose = false;
end

if ~exist('template','var')
    useTemplate = false;
else
    useTemplate = true;
end

if ndims(imageIn) == 4
    [yPixels, xPixels, nWindows, nTrials] = size(imageIn);
    Y = reshape(imageIn,[yPixels, xPixels, nWindows*nTrials]);
elseif ndims(imageIn) == 3
    [yPixels, xPixels, nWindows] = size(imageIn);
    Y = imageIn;
end

Y = single(Y);                 % convert to single precision
Y = Y - min(Y(:));

%% set parameters (first try out rigid motion correction)
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200, 'shifts_method', 'cubic');
options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200, 'shifts_method', 'cubic');

%% perform motion correction
if useTemplate
    tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid,template); toc %#ok<*ASGLU>
else
    tic; [M1,shifts1,template1,options_nonrigid] = normcorre_batch(Y,options_nonrigid); toc
end
%% non-rigid motion correction (or whatever comparison you want)
if verbose
    disp('also testing non-rigid motion correction (verbose is selected)')
    tic; [M2,shifts2,template2,options_rigid] = normcorre(Y,options_rigid); toc %#ok<*ASGLU>
end

%% compute metrics & plot the shifts/video of the result
if verbose
    nnY = quantile(Y(:),0.005);
    mmY = quantile(Y(:),0.995);
    
    [cY,mY,vY] = motion_metrics(Y,10);
    [cM1,mM1,vM1] = motion_metrics(M1,10);
    [cM2,mM2,vM2] = motion_metrics(M2,10);
    T = length(cY);
    
    % values for prettier plotting
    mYp = normalizeFrame(mY, min(mY(:)), max(mY(:)));
    mM1p = normalizeFrame(mM1, min(mM1(:)), max(mM1(:)));
    mM2p = normalizeFrame(mM2, min(mM2(:)), max(mM2(:)));
    
    % plot metrics
    figure;
    ax1 = subplot(2,3,1); imagesc(mYp,[0,1]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,3,2); imagesc(mM2p,[0,1]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    ax3 = subplot(2,3,3); imagesc(mM1p,[0,1]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','non-rigid','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
    xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
    xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2,ax3],'xy')
    % plot shifts
    
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
    shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
    shifts_x = squeeze(shifts_nr(:,1,:))';
    shifts_y = squeeze(shifts_nr(:,2,:))';
    
    patch_id = 1:size(shifts_x,2);
    str = strtrim(cellstr(int2str(patch_id.')));
    str = cellfun(@(x) ['patch # ',x],str,'un',0);
    
    figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','non-rigid','rigid');
    title('correlation coefficients','fontsize',14,'fontweight','bold')
    set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2);
    title('displacements along x','fontsize',14,'fontweight','bold')
    set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2);
    title('displacements along y','fontsize',14,'fontweight','bold')
    xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')
    
    
    %% plot a movie with the results
    saveBool = input('Would you like to create a video? (y/n): ','s');
    if strcmp(saveBool,'y')
        h1 = figure();
        set(h1, 'Position', [50 50 1800 700])
        minAct = min(Y(:));
        maxAct = max(Y(:));
        minActM2 = min(M1(:));
        maxActM2 = max(M1(:));
        frameNum = 1;
        for t = 1:1:T
            rawFrame = normalizeFrame(Y(:,:,t), minAct, maxAct);
            correctedFrame = normalizeFrame(M1(:,:,t), minActM2, maxActM2);
            
            subplot(121);imagesc(rawFrame,[0 1]);
            xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
            title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
            
            subplot(122);imagesc(correctedFrame,[0 1]);
            xlabel('rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
            title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
            
            set(gca,'XTick',[],'YTick',[]);
            drawnow;
            pause(0.01);
            
            %if ~mod(t,30) % save every 30th frame for displaying to video
            frames(frameNum) = getframe(gcf);   %#ok<AGROW>
            frameNum = frameNum+1;
            %end
        end
        
        %% save video
        startDir = pwd;
        saveDir = uigetdir();
        cd(saveDir);
        fprintf('Saving...')
        
        vidObj = VideoWriter(strcat('normcore_registration',date),'MPEG-4');
        vidObj.FrameRate = 30;
        open(vidObj);
        writeVideo(vidObj,frames);
        close(vidObj);
        
        cd(startDir)
        fprintf('Done.\n')
    end
end

%% pass back the result
disp('returning rigid-body motion corrected data')
if ndims(imageIn) == 4
    imageOut = reshape(M1,[yPixels, xPixels, nWindows, nTrials]);
else
    imageOut = M1;
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%% sub functions %%%%%%%%%%%%%%%%%%%%%%%%%%%

function frameOut = normalizeFrame(frameIn, minAct, maxAct)
% normalise data [0 1]
frameOut = (frameIn - minAct)./(maxAct-minAct);
% nonlinear scale between [0 1]
frameOut = nthroot(frameOut,4);
end
