function [B,L,CC,stats] = getROIprocess(img,roiSz)
DbFlag = 0;
% roiSz = 20;
DOGsz = [0.1 1]*roiSz;
SizeOfSmallParticleToRemove = round(pi*(roiSz/2*0.75)^2);
SizeOfLargeParticleToRemove = round(pi*(roiSz*0.75)^2);
% t = A(:,:,1);
% figure,imagesc(img);
% ----- Apply DOG filter to enhance activated pixels ----------------------
tmf = medfilt2(img);
tf = imgaussfilt(tmf,DOGsz(1));
if DbFlag,figure(1000000+10),imshow(tf,[ ],'colormap',parula(256));end
tfbg = imgaussfilt(img,DOGsz(2));
if DbFlag,figure(1000000+20),imshow(tfbg,[0,1],'colormap',parula(256));end
td = tf-tfbg;
if DbFlag,figure(1000000+30),imshow(td,[0 0.5],'colormap',parula(256));end

% ----Triangle threshold---------------------------------------------------
%     figure(40);
%     h = histogram(td,50);
%     hold on;
[N,edges] = histcounts(td,50);
[Vmax,Imax] = max(N);
Xmax = mean(edges(Imax:(Imax+1)));
Vtail = N(end);
Xtail = mean(edges((end-1):end));
%plot([Xmax Xtail],[Vmax,Vtail]);
% point2line distance
dA = zeros(length(Imax+1:length(N)),1);
ptxA = dA;
for pti = 1:length(dA)
    ptx = mean(edges((Imax:(Imax+1))+pti));
    pty = N(Imax+pti);
    va = [Xmax-Xtail,Vmax-Vtail,0];
    vb = [ptx-Xtail,pty-Vtail,0];
    c = cross(va,vb);
    d = norm(c)./norm(va);
    dA(pti) = d;
    ptxA(pti) = ptx;
end
% figure,plot(dA);
[~,Imax2] = max(dA);
th = ptxA(Imax2);
%     plot(th,Values(Imax+Imax2),'r*');
%     hold off
%--------- binary map and extract ROIs-------------------------------------
mask = td>th;
maskc0 = bwareaopen(mask,SizeOfSmallParticleToRemove);%remove small particals
maskc1 = bwareaopen(mask,SizeOfLargeParticleToRemove);%keep large particals
maskc = maskc0-maskc1;%remove larege paritcals
if DbFlag,figure(1000000+50),imshow(maskc,[0 1],'colormap',parula(256));end
maskc = imfill(maskc,'holes');
[B,L,~,~] = bwboundaries(maskc);
if DbFlag, figure(1000000+60),imshow(L,[],'colormap',parula(256)); hold on; end
% --------- remove non-round objects ---------------------------------------
stats = regionprops(L,'Area','Centroid');
threshold = 0.5;
for k = 1:length(B)
  boundary = B{k};
  delta_sq = diff(boundary).^2;    
  perimeter = sum(sqrt(sum(delta_sq,2)));
  area = stats(k).Area;
  metric = 4*pi*area/perimeter^2;
  if metric<threshold
      L(L==k)=0;
  end
  if DbFlag,text(stats(k).Centroid(1),stats(k).Centroid(2),num2str(metric));end
end
% ---------- re-trace boundaries ------------------------------------------
[B,L,~,~] = bwboundaries(logical(L));
if DbFlag, figure(1000000+70),imshow(L,[],'colormap',parula(256)); end
% if DbFlag, pause(0.5);end
% ----------- ROI pixel list ---------------------------------------------
CC = bwconncomp(logical(L));
stats = regionprops(L,'Area','Centroid');