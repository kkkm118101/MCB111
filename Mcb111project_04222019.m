MCB111_project

%% Morphing
%Resize the image
RefImage=imresize(imread('BEAM_Emx1Cre_WT_Ctip1_hom_2_Slide_6.tif'), [786 1048]);
Tomorph=imresize(imread('01182019_BEAM_Ctip1_null_E15_P7_mFP_4_Slide_4.tif'), [786 1048]);

subplot(1,2,1); imshow(RefImage);title('RefImage')
subplot(1,2,2); imshow(Tomorph); title('Tomorph')
imwrite(RefImage, 'RefImage.tif');
imwrite(Tomorph, 'Tomorph.tif');

%Determine all points for morphing
captureWarpPoints('BEAM_Emx1Cre_WT_Ctip1_hom_2_Slide_6.tif');
captureWarpPoints('01182019_BEAM_Ctip1_null_E15_P7_mFP_4_Slide_4.tif');
    % Capture all points (default: 25 points)
    % Double-click and the points will be saved in a file named "filename".pts.mat 
    % Repeat this process with the second image (tomorph.tif in our example).

%Morphing
warpRunner('Tomorph.tif','01182019_BEAM_Ctip1_null_E15_P7_mFP_4_Slide_4.tif.pts.mat', 'BEAM_Emx1Cre_WT_Ctip1_hom_2_Slide_6.tif.pts.mat','tif'); %Warp tomorph into reference
    % Do morphing for each channel by changing line 13 of warpImage.m, from (:,:,1) to (:,:,2) and
    % line 30 of warpRunner.m, from 'warped1_' to 'warped2_'
warpRunner('Tomorph.tif','01182019_BEAM_Ctip1_null_E15_P7_mFP_4_Slide_4.tif.pts.mat', 'BEAM_Emx1Cre_WT_Ctip1_hom_2_Slide_6.tif.pts.mat','tif'); %Warp tomorph into reference

subplot(1,2,1);imshow('warped1_Tomorph.tif');title('warped1_Tomorph.tif');
subplot(1,2,2);imshow('warped2_Tomorph.tif');title('warped2_Tomorph.tif');

%% Save previously adjusted image with Image viewer and export them in worksplace as adj_Ch1_wTomorph, Ch2_wTomorph...
imwrite(Ch1_wTomorph, 'adj_Ch1_wTomorph.tif');
imwrite(Ch2_wTomorph, 'adj_Ch2_wTomorph.tif');
%% Select region of interest and extract it

%Channel 1
%Crop image
grayImage = imread('adj_Ch1_wTomorph.tif');
imshow(grayImage);
h = imrect;
position = wait(h);
croppedImage = imcrop(grayImage, position);
figure;
imshow(croppedImage);
imwrite(croppedImage,'ROI_adCh1_wTomorph.tif','tif');

%Selection and extraction of ROI
[X,MAP]=imread('ROI_adCh1_wTomorph.tif');
imshow('ROI_adCh1_wTomorph.tif');

lsz=length(size(X));
if lsz==3
    X=rgb2gray(X); 
end

im0=X;
X=double(X);
h_im=imshow(im0);title('Choose ROI from the figure');
uiwait(msgbox('Draw a polygon. double-click to accept/finish it'))
h=impoly;
polyPoints=h.getPosition;
position=wait(h);
e = impoly(gca,position);
BW=createMask(h,h_im);
ROI=X.*BW;
figure, imshow(ROI, []);
im_DIF=X-ROI;
X=X-ROI;
figure, imshow(im_DIF, []);
ROI=uint8(ROI); im_DIF=uint8(im_DIF); 
imwrite(ROI, 'ROIex_adCh01_wTomorph.tif', 'tif');
imwrite(im_DIF,'Image difference_ROIex_adCh01_wTomorph.tif', 'tif');

%% Delauney triangulation
% Copy coordinates of ROI position by clicking right button of mouse and
% put values in a vector 'Coord'.
Coord=[28.0000000000001 248;20.0000000000001 189;23.0000000000001 77.0000000000001;58 39.0000000000001;127 28.0000000000001;222 33.0000000000001;305 68.0000000000001;378 150;410 183;466 212;481 285;483 392;468 485;414 575;343 533;394 446;416 395;418 342;411 303;399 284;360 254;322 234;253 202;233 191;210 195;188 198;157 216;93 247]

%Delaunay triangulation
T=delaunayTriangulation(Coord(:,1), Coord(:,2));
hold on; triplot(T);
%Compute the in/out status.
I0=T.isInterior;
patch('faces',T(I0,:), 'vertices',T.Points, 'FaceColor','c');
axis equal;
% Cannot eliminate triangles outside of ROI. Possibly use convex hull 2D
% Put into rectangles
tform=[0 0 0 0 0 0 0 0 0 0 0 0 .25 .5 .75 1 1.25 1.5 1.75 8 8 8 8 8 1.5 1 .5 0 ; 
       1 0.9 0.8 0.75 0.7 0.6 0.5 0.4 0.3 0.7 .25 0 0  0 0 0 0 0 0 0 .25 .5 .75 1 1 1 1 1];
tform=tform';
warpRunner('ROIex_adCh01_wTomorph.tif',Coord,tform,'tif')
t=fitgeotrans('ROIex_adCh01_wTomorph.tif',tform,'pwl');

%% Channel 2
%Crop image
grayImage = imread('adj_Ch2_wTomorph.tif');
imshow(grayImage);
h = imrect;
position = wait(h);
croppedImage = imcrop(grayImage, position);
figure;
imshow(croppedImage);
imwrite(croppedImage,'ROI_adCh2_wTomorph.tif','tif');

%Selection and extraction of ROI
[X,MAP]=imread('ROI_adCh2_wTomorph.tif');
imshow('ROI_adCh2_wTomorph.tif')
lsz=length(size(X));
if lsz==3
    X=rgb2gray(X); 
end

im0=X;
X=double(X);
h_im=imshow(im0);title('Choose ROI from the figure');
uiwait(msgbox('Draw a polygon. double-click to accept/finish it'))
h=impoly;
polyPoints=h.getPosition;
position=wait(h);
e = impoly(gca,position);
BW=createMask(h,h_im);
ROI=X.*BW;
figure, imshow(ROI, []);
im_DIF=X-ROI;
X=X-ROI;
figure, imshow(im_DIF, []);
ROI=uint8(ROI); im_DIF=uint8(im_DIF); 
imwrite(ROI, 'ROIex_adCh02_wTomorph.tif', 'tif');
imwrite(im_DIF,'Image difference_ROIex_adCh02_wTomorph.tif', 'tif');

%% Delauney triangulation
% Copy coordinates of position by clicking right button of mouse and
% make a vector 'Coord'.
Coord=[28.0000000000001 248;20.0000000000001 189;23.0000000000001 77.0000000000001;58 39.0000000000001;127 28.0000000000001;222 33.0000000000001;305 68.0000000000001;378 150;410 183;466 212;481 285;483 392;468 485;414 575;343 533;394 446;416 395;418 342;411 303;399 284;360 254;322 234;253 202;233 191;210 195;188 198;157 216;93 247]

% Delaunay triangulation
T=delaunayTriangulation(Coord(:,1), Coord(:,2));
hold on; triplot(T);
%Compute the in/out status.
I0=T.isInterior;
patch('faces',T(I0,:), 'vertices',T.Points, 'FaceColor','c');
axis equal;

%Put into rectangle
tform=[0 0 0 0 0 0 0 0 0 0 0 0 .25 .5 .75 1 1.25 1.5 1.75 8 8 8 8 8 1.5 1 .5 0 ; 
       1 0.9 0.8 0.75 0.7 0.6 0.5 0.4 0.3 0.7 .25 0 0  0 0 0 0 0 0 0 .25 .5 .75 1 1 1 1 1];
tform=tform';
strCoord=string(Coord);
strtform=string(tform);
warpRunner('ROIex_adCh01_wTomorph.tif',Coord,tform,'tif')
t=fitgeotrans('ROIex_adCh01_wTomorph.tif',tform,'pwl');
%final=imwarp(RefImage_1, t);
%figure(),imshow(final);

%% Binarize it
BWCh2_wTomorph= imbinarize(imread('adj_Ch2_wTomorph.tif'),'adaptive','Sensitivity', 0.55);
imshowpair(imread('adj_Ch2_wTomorph.tif'),BWCh2_wTomorph,'montage')
imwrite(BWCh2_wTomorph,'BWCh2adj_wTomorph.tif');
imshow('BWCh2adj_wTomorph.tif')


%% Rotate
I1=imread('croppedROI_Ch1_tdt_Tomorph.tif');
J1 = imrotate(I1,45);imshow(J1);title('J1')
imwrite(J1,'RcroppedROI_Ch1_tdt_Tomorph.tif');
I2=imread('croppedROI_Ch2_gfp_Tomorph.tif');
J2 = imrotate(I2,45); imshow(J2);title('J2')
imwrite(J2, 'RcroppedROI_Ch2_gfp_Tomorph.tif');

%% Measure intensity
nbins=400;  
binrange= linspace(0,400,nbins);   %Allow for 16 bit, but round up so each bin is 66 intensity levels.
[thisCounts1,ind] = histc(J1/sum(sum(J1)), binrange); % Get counts for this image.
figure(), subplot(2,2,1); plot(binrange, thisCounts1, 'k','linewidth', 2);title('Contralateral distribution of Ctip1 null neurons');xlabel('bins');ylabel('intensity');
subplot(2,2,2); imshow('RcroppedROI_Ch1_tdt_Tomorph.tif'); title('Ctip1 null neurons projection in contralateral cortex');
[thisCounts2,ind] = histc(J2/sum(sum(J2)), binrange); % Get counts for this image.
subplot(2,2,3); plot(binrange,thisCounts2, 'k', 'linewidth', 2);title('Contralateral distribution of wt neurons');xlabel('bins');ylabel('intensity');
subplot(2,2,4); imshow('RcroppedROI_Ch2_gfp_Tomorph.tif'); title('WT neurons projection in contralateral cortex')

%% Visualize contralateral hemispheric projection.
figure(); subplot(1,2,1); contour(J1);title('Contralateral distribution of Ctip1 null neurons');
subplot(1,2,2); contour(J2);title('WT neurons projection in contralateral cortex')