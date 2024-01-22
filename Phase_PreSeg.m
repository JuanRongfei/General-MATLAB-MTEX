%% calculate grains based on raw data
[grains_raw,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',15*degree);
figure;
plot(grains_raw, 'coordinates','on');
hold off;
saveas(gcf, 'grains_raw.png')
% Raw GOS map plot
figure;
plot(grains_raw,grains_raw.GOS./degree,'coordinates','on')
% mtexColorbar('Location','westoutside','FontSize',24,'FontWeight','bold');
% caxis([0,10])
hold off;
saveas(gcf, 'grains_raw_GOS.png')

% smooth the data by using the median Filter to smoothe the data
F = medianFilter;
% define the size of the window to be used for finding the median
F.numNeighbours = 1;
ebsd_smoothed = smooth(ebsd('indexed'),F,'fill');
ebsd_s= ebsd_smoothed;
clear ebsd_smoothed
ebsd_smoothed= ebsd_s;

% Plot smoothed the EBSD data 
figure;
plot(ebsd_smoothed);
saveas(gcf, 'ebsd_smoothed1.png')
% Plot smoothed the EBSD data，damaged/voids color in black
figure;
plot(ebsd_smoothed('notIndexed'),'FaceColor','black');
hold on;
plot(ebsd_smoothed);
hold off;
saveas(gcf, 'ebsd_smoothed2.png')
% Plot KAM
figure;plot(ebsd_smoothed,ebsd_smoothed.KAM('threshold',5*degree,'order',1)./degree,'micronbar','on');
figure;plot(ebsd,ebsd.KAM('threshold',5*degree,'order',1)./degree,'micronbar','on');
setColorRange([0 3],'current');
mtexColorbar location eastoutside
mtexColorMap jet
hold on
plot(grains_Bcc.boundary,'lineWidth',1,'lineColor','w')
hold off
saveas(gcf, 'EBSD2cRange_KAM5MapO1+GB15_allw.tif')

%% calculating the grains variable based on the smoothed data
[grains_smoothed,ebsd_smoothed.grainId,ebsd_smoothed.mis2mean] = calcGrains(ebsd_smoothed,'angle',15*degree);
figure;
plot(grains_smoothed, 'coordinates','on');
saveas(gcf, 'grains_smoothed1.png')
% Plot smoothed grains，damaged/voids color in black
figure;
plot(grains_smoothed('notIndexed'),'FaceColor','black');
hold on;
plot(grains_smoothed, 'coordinates','on');
hold off;
saveas(gcf, 'grains_smoothed2.png')

% remove some very small grains, containing less than 5 pixel grains
ebsd_smoothed(grains_smoothed(grains_smoothed.grainSize<5)) = [];
[grains_5,ebsd_s.grainId,ebsd_s.mis2mean] = calcGrains(ebsd_s,'angle',15*degree);
% smooth grain boundaries
grains_5 = smooth(grains_5,2);

% plot smooth grains with boundaries.
figure;
plot(grains_5,'coordinates','on')
saveas(gcf, 'grains_5.png')
figure;
plot(grains_5,grains_5.GOS./degree,'coordinates','on')
mtexColorbar('Location','westoutside','FontSize',24,'FontWeight','bold');
caxis([0,10])
saveas(gcf, ['grains_5_GOS.png'])

%% Calculate the Grain average IQ value
threshold = 5;
% gridify the raw ebsd data 
ebsd_grid = ebsd.gridify;
% calculating the grains
[grains,ebsd_grid.grainId,ebsd_grid.mis2mean] = calcGrains(ebsd_grid,'angle',15*degree);
% define grains with grain size larger than 5 pixels 
ebsd_condition = ebsd_grid(grains(grains.grainSize>5));
% calculating the grains containing more than 5 pixels
[grains_5pixel,ebsd_condition.grainId,ebsd_condition.mis2mean] = calcGrains(ebsd_condition,'angle',15*degree);
% calculating the
geo=ebsd_condition.grainId;
% iq_mean only have the data from ebsd_condition, grains_5pixel
iq_mean=zeros(length(geo),1);
kam_mean=zeros(length(geo),1);
bs_mean=zeros(length(geo),1);

for i=1:max(grains_5pixel.id)
    meanIQ=mean(ebsd_condition(grains_5pixel(i)).bc);
    meanKAM=mean(ebsd_condition(grains_5pixel(i)).KAM);
    meanBS=mean(ebsd_condition(grains_5pixel(i)).bs);
    k=find(geo==i);
        for j=1:length(k)
        iq_mean(k(j))=meanIQ;
        kam_mean(k(j))=meanKAM;
        bs_mean(k(j))=meanBS;
        end
end

y=ebsd_grid.y; y=y(:);
x=ebsd_grid.x; x=x(:);
A_1=[x,y];
x2=ebsd_condition.x;
y2=ebsd_condition.y;
B=[x2,y2,iq_mean];
ebsd_condition_2 = ebsd_grid(grains(grains.grainSize<=5));
x3=ebsd_condition_2.x;
y3=ebsd_condition_2.y;
iq_condition_2 = ebsd_condition_2.bc;
kam_condition_2 = ebsd_condition_2.KAM;
C1=[x3,y3,iq_condition_2];
C2=[x3,y3,kam_condition_2];
D1=[B;C1];
D2=[B;C2];
E1=sortrows(D1,[1,2]);
E2=sortrows(D2,[1,2]);
% Create the column vector iq_final from the elements of row 3 of matrix a: iq_final=E(:,3)
iq_final=E1(:,3);[m,n]=size(ebsd_grid.x);
iq_final=reshape(iq_final,[round(m),round(n)]);
ebsd_grid.prop.GrainAverIQ = iq_final; 
kam_final=E2(:,3);[m,n]=size(ebsd_grid.x);
kam_final=reshape(kam_final,[round(m),round(n)]);
ebsd_grid.prop.GrainAverKAM = kam_final;
deltaX=ebsd_grid.x(1,2)-ebsd_grid.x(1,1);
dx=ebsd_grid.x(end)-ebsd_grid.x(1)+ deltaX;
y=ebsd_grid.y; y=y(:);
x=ebsd_grid.x; x=x(:);
deltaY=y(2,1)-y(1,1);
dy=ebsd_grid.y(end)-ebsd_grid.y(1)+deltaY;
ind = inpolygon(ebsd_grid,[ebsd_grid.x(1),ebsd_grid.y(1),dx,dy]);
ebsd = ebsd_grid(ind);
plot(ebsd,ebsd.GrainAverIQ,'figSize','large')
mtexColorMap
hold on;
plot(ebsd,ebsd.GrainAverKAM,'figSize','large','faceAlpha',0.7);
hold off;
%calculating the size of the ebsd data
ebsd_sizeY=dy/deltaY;
ebsd_sizeX=dx/deltaX;
ind_1 = ebsd ('Iron bcc (old)').id;
ind_2 = ebsd ('notIndexed').id;
ind_bcc_NI = [ind_1; ind_2];
ebsd_bcc_NI = ebsd (ind_bcc_NI);

%% calculation GOS of ebsd 
Gos = zeros(size(ebsd_s.grainId));
% Create an all-0 matrix, but with the same shape as the ebsd
for i=1:length(ebsd_s.grainId)
   
    gos_mean = grains_smoothed(ebsd_s.grainId(i)).GOS; 
   
    Gos(i) = gos_mean;
end

% The first line defines how many steps the loop takes and the loop condition
% The second line of the loop body, what the loop is going to do
% assign the calculated Gos values to ebsd_wo_M properties
ebsd_s.prop.GOS = Gos;
ebsd.prop.GOS = Gos;
ebsd_bcc = smooth(ebsd_bcc_NI('indexed'),F,'fill');
ebsd_bcc.prop.GOS = Gos;


%% plot threshold
BChistBcc = ebsd_bcc.GrainAverIQ;
figure;
ksdensity(BChistBcc);
saveas(gcf, 'BChistBcc.png')
GOShistBcc = ebsd_bcc.GOS./degree;
figure;
ksdensity(GOShistBcc);
saveas(gcf, 'GOShistBcc.png')