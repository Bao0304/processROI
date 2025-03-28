%% load data
fnA=dir('*.tif');
dataLabel = 'astrocyte';
roiSz = 4;
botA = cell(size(fnA,1),1);
for mi = 1:size(fnA,1)
    fn = fnA(mi).name;
    disp(fn);
    tiff_info = imfinfo(fn);
    frameNum = length(tiff_info);
    if frameNum>120
        frameNum=150;
    end
    A = zeros(tiff_info(1).Height,tiff_info(1).Width,frameNum);
    for fi = 1:frameNum
        A(:,:,fi) = imread(fn,fi);
    end
    refImg = mean(A,3);
    
        %% find max lum 5 data points of a pixel.
    B = sort(A,3,'descend');
    Bm = mean(B(:,:,1:10),3);
  
    %% get ROI
    LA = zeros(size(A));
    BA = cell(frameNum,1);
    CCA = cell(frameNum,1);
    statsA = cell(frameNum,1);
    for fi = 1:frameNum
        [B,L,CC,stats] = getROIprocess(A(:,:,fi),roiSz);
        BA{fi} = B;
        LA(:,:,fi) = L;
        CCA{fi} = CC;
        statsA{fi} = stats;
        fprintf('frame %d, number of roi %d.\n', fi,max(L(:)));
    end
    %% conbine data, and remove overlap cells
    clear CC
    clear B
    clear stats
    B = BA{1};
    CC = CCA{1}.PixelIdxList;
    stats = statsA{1};
    for fi = 2:frameNum
        B = cat(1,B,BA{fi});
        CC = cat(2,CC,CCA{fi}.PixelIdxList);
        stats = cat(1,stats,statsA{fi});
    end
    centroids =  extractfield(stats,'Centroid');
    area =  extractfield(stats,'Area');%area>400 will be discarded later
    centroids = reshape(centroids,2,[]); %xy 2xn
    [cellpos,T]=MergeClosePoints(centroids,roiSz*1.5);
    %% merge pixel list of CC
    SizeOfSmallParticleToRemove = round(pi*(roiSz/2*0.9)^2);
    CCB = {};
    cci =1;
    for ci = 1:max(T(:))
        mz = zeros(size(Bm));
        cellids=(find(T==ci));
        pl = cat(1,CC{cellids});
        for pli = 1:size(pl,1)
            mz(pl(pli)) = mz(pl(pli)) + 1;
        end
        M = max(mz(:));
        pth = floor(M*0.75);
        
        mz(mz<=pth)=0;
                % -------------- remove small rois ----------------
        pn = sum(logical(mz(:)));
        if pn<SizeOfSmallParticleToRemove
            continue;%CCB{ci}=[];
        end
        % -------------- remove non round rois ------------
        [B,L] = bwboundaries(logical(mz));
        stats = regionprops(L,'Area','Centroid');
        threshold = 0.5;
        boundary = B{1};
        delta_sq = diff(boundary).^2;
        perimeter = sum(sqrt(sum(delta_sq,2)));
        area = stats(1).Area;
        metric = 4*pi*area/perimeter^2;
        if metric<threshold
            continue;%CCB{ci}=[];
        end
        CCB{cci}=find(mz>pth);
        cci = cci + 1;
    end
    roiNum = length(CCB);
    
    Ln = zeros(size(Bm));
    for cci=1:roiNum
        p = CCB{cci};
        Ln(p) = cci;
        [y,x] = ind2sub(size(Bm),p);
        plot(x,y,'.');
        %     disp(length(x));
    end
   
    %% BOT
    roiNum = length(CCB);
    bot = zeros(roiNum,size(A,3));
    Ar = reshape(A,size(A,1)*size(A,2),[]);
    for cci = 1:roiNum
        p = CCB{cci};
        bot(cci,:) = mean(Ar(p,:));
    end
    
    botA{mi} = bot;
end
save botA_ROI botA
%% combine data into single matrix, 
botB = [];
for mi = 1:size(fnA,1)
    bot = botA{mi};
    %if size(bot,2)==120%in case of 120, discard first 40
       % bot(:,1:40)=[];
   % end
    %if size(bot,2)==80% in case of 80, binning to 40.
   %     botr = reshape(bot,size(bot,1),2,40);
   %     bot = squeeze(mean(botr,2));
   % end
    botB = cat(1,botB,bot);
end
save botB_ROI botB
totalRoiNum = size(botB,1);
%%  show botB
load botB_ROI
t = 1:40;

imagesc(botB);
yticks(1:10:totalRoiNum);
ylabel('roi ID');
xlabel('Time');
%% botB Normalize
botN = zeros(size(botB));
for ri = 1:totalRoiNum
%     botZ(ri,:) = zscore(botB(ri,:));
    tt = botB(ri,:);
    ttn = (tt-min(tt))/(max(tt)-min(tt));
    botN(ri,:) = ttn;
end

%%

%% sort normalized botB, show it.
% Z = linkage(botZ,'ward');
% [~,~,I]=dendrogram(Z,0);%'ColorThreshold',cutoff
botNm = mean(botN(:,11:40),2);
[~,I] = sort(botNm,'descend');
botC = botN(I,:);
cmap = gray(256);
cmap = flipud(cmap);
figure,imagesc(botC); colormap(cmap);
xlabel('Time (s)');
ylabel('ROI ID');
title('Normalized traces of ROIs');
saveas(gcf,['RasterPlot_' dataLabel '_ROI'],'epsc');
saveas(gcf,['RasterPlot_' dataLabel '_ROI'],'bmp');

%% raster plot
%% bar plot

botNm = mean(botN);
% botbmr = reshape(botBm,2,[]);
% m = mean(botbmr,1);
m = mean(botNm,1);
bar(m,'k');
xlabel('Time (s)');
ylabel('Average of normalized roi traces');
saveas(gcf,['BarPlot_' dataLabel '_ROI'],'epsc');
saveas(gcf,['BarPlot_' dataLabel '_ROI'],'bmp');

save botNm_ROI botNm
%%
% epscFn = 'test';
% saveas(gcf,epscFn,'epsc2');