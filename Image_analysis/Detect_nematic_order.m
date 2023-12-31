%% Image_analysi_algorithm to detect the nematic director
% Copyright (C) <2023>  <Francine Kolley>
% The code was written for multi-channel images of insect flight muscle
% Pipeline in collaboration with Benoit Dehapiot
% Francine Kolley
% Physics of Life, Benjamin M. Friedrich group
% TU_dresden 
% contact: francine.kolley@tu-dresden.de
% Latest code 07-2022

% STEPS:
% (1) run this script
% (2) rund Correlation_function.m

clear all
close all

% #########################################################################
%                             PARAMETER  
% #########################################################################

% Background subtraction
GBlur_Sigma = 2; 

% Mask treshold to filter interesting structures for ROI 
Mask_Thresh =350; %it varies for different time points of myofibrillogenesis

% ROIsMask (obtained from Mask)
ROI_Size = 30; %[pixels] 
ROI_Thresh = 0.25; % minimal value used, some stages increase up to 0.85

% Tophat filtering
Tophat_Sigma = 28; 

% Steerable Filter
Steerable_Sigma = 2; % Size for Steerable filter  

% Padding
Pad_Angle = 6; % Window size --> we will crop image to ensure linescans are not reaching edge of image
Pad_Corr = 3; % [pixel]
ROI_Pad_Corr = (ROI_Size+(ROI_Size*(Pad_Corr-1))*2);


% #########################################################################
%                             IMAGE READ IN 
% #########################################################################

input_figures_directory='insert_path\';
image_filename=['insert_image_name.tif'];

% define the order of the different channels, number of channels
[actinin_image,actin_image,myosin_image,titin_image]=read_in_images(sprintf('%s\\%s',input_figures_directory,image_filename),4);

% save the path for meta_data later 
RawPath = strcat(input_figures_directory,'\',image_filename);
% Open Data and read in information
info = imfinfo(RawPath);
% define the slice you want to observe from the 3D image
slice=1;
Raw =titin_image(:,:,slice);

% Get the image_factor 
factor=info.XResolution;
pixSize=1/factor; %factor to go from [pixel] from [\mu m] 
clear info 

% Get variables
nY = size(Raw,1);
nX = size(Raw,2);  
qLow = quantile(Raw,0.001,'all');
qHigh = quantile(Raw,0.999,'all');

% #########################################################################
%                 MASK AND BACKGROUND SUBSTRACTION  
% #########################################################################
choice = 1;
while choice > 0   

    % Crop Raw image
    nGridY = floor(nY/ROI_Size);
    nGridX = floor(nX/ROI_Size);
    Raw(nGridY*ROI_Size+1:end,:) = [];
    Raw(:,nGridX*ROI_Size+1:end) = [];
    nYCrop = size(Raw,1);
    nXCrop = size(Raw,2);
       
    % Substract background
    Raw_BGSub = imgaussfilt(Raw,GBlur_Sigma); % Gaussian blur #2
    

    % Create Mask
    % Binary mask; I<Threshold --> Black, I>Threshold --> White
    Raw_Mask = Raw_BGSub;
    Raw_Mask(Raw_Mask<Mask_Thresh) = 0;
    Raw_Mask(Raw_Mask>=Mask_Thresh) = 1;

    % Create ROIsMask
    Raw_ROIsMask = zeros(nGridY,nGridX);
    
    for i=1:nGridY
        for j=1:nGridX
            temp = mean(Raw_Mask(ROI_Size*i-(ROI_Size-1):ROI_Size*i,ROI_Size*j-(ROI_Size-1):ROI_Size*j));
             white_percentage(i,j)=mean(temp);%how many pixel per ROI white
            if mean(temp(:)) > ROI_Thresh 
                if i >= Pad_Angle && i <= nGridY-(Pad_Angle-1) && j >= Pad_Angle && j <= nGridX-(Pad_Angle-1)
                    Raw_ROIsMask(i,j) = 1;
                   
                end

            end
        end
    end
    
% Display .................................................................
    
    % Mask
    subplot(2,1,1) 
    imshow(Raw_Mask,[0 1])
    title(strcat(...
        'Mask (thresh. =',{' '},num2str(Mask_Thresh),')'));
    
    % ROIsMask
    subplot(2,1,2) 
    imshow(Raw_ROIsMask,[0 1])
    title(strcat(...
        'ROIsMask (ROIs size = ',{' '},num2str(ROI_Size),{' '},'pix. ;',...
        {' '},'ROIs thresh. =',{' '},num2str(ROI_Thresh),{' '},';'));
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %pause(10)
% Dialog box ..............................................................
    
    choice = questdlg('What next?', ...
        'Menu', ...
        'Modify Parameters','Proceed','Proceed');
    switch choice
        case 'Modify Parameters'
            choice = 1;
        case 'Proceed'
            choice = 0;
    end

    if choice == 1

        prompt = {
            'Mask_Thresh :',...
            'ROI_Thresh :'};

        definput = {
            num2str(Mask_Thresh),...
            num2str(ROI_Thresh)
            };
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
 
        Mask_Thresh = answer(1,1); 
        ROI_Thresh = answer(2,1); 
        
        close
    end

    if choice == 0
        
        close
    end
    
end

% #########################################################################
%                         TOP HAT TRANSFORMATION  
% #########################################################################
% top hat transormation with fixed values for sigma! 
Raw_Tophat = imtophat(Raw_BGSub,strel('disk',Tophat_Sigma)); 
    

% #########################################################################
%                   MEASURE LOCAL NEMATIC ORDER
% #########################################################################
choice = 1;
while choice > 0   

    % Steerable Filter
    RawRSize = double(imresize(Raw,[nGridY nGridX],'nearest'));
    MaskRSize = double(imresize(Raw_Mask,[nGridY nGridX],'nearest'));
    [res,~,~,rot] = steerableDetector(RawRSize,2,Steerable_Sigma,180);
    for i=1:size(rot,3)
        temp = rot(:,:,i);
        temp(MaskRSize==0) = NaN;
        rot(:,:,i) = temp;
    end

    % Make AngleMap
    AngleMap = zeros(nGridY,nGridX);
    parfor i=1:nGridY
        for j=1:nGridX
            if Raw_ROIsMask(i,j) == 1            
                Crop = rot(i-(Pad_Angle-1):i+(Pad_Angle-1),j-(Pad_Angle-1):j+(Pad_Angle-1),:);
                idxMax = NaN(size(Crop,3),1);
                for k=1:size(Crop,3)
                    temp = Crop(:,:,k);
                    idxMax(k,1) = nanmean(temp(:));
                end
                [M,I] = max(idxMax);
                AngleMap(i,j) = I;            
            end
        end
    end



% Display .................................................................
    % Display angles
  

     % Raw for control
     subplot(3,1,1)
     imshow(Raw,[qLow qHigh])
     title('Raw')

    % Steerable filter
    subplot(3,1,2)
    imshow(res,[min(res(:)) max(res(:))])
    title(strcat(...
        'Steerable filter (sigma =',{' '},num2str(Steerable_Sigma),{' '},'pix.)'));


     %nematic
      Raw_new=ones(size(Raw,1),size(Raw,2)); %just to get a black background
      subplot(3,1,3); hold on
      imshow(Raw_new,[qLow qHigh])
      hold on
for i=1:nGridY
    for j=1:nGridX
        if Raw_ROIsMask(i,j) == 1  % if ROi is white = of interest
        
            Angle = AngleMap(i,j); % find right angle for window
           
                k=0;
                xi = ROI_Size*j+ROI_Size *[-1 1]*cos((90-Angle)*-1*pi/180)+k*cos((Angle)*pi/180);
                yi= ROI_Size*i+ROI_Size *[-1 1]*sin((90-Angle)*-1*pi/180)+k*sin((Angle)*pi/180);
            
            plot(xi,yi,'c')
            hold on
               
            else
        end
    end
end
    

    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
   
% Dialog box ..............................................................
    
    choice = questdlg('What next?', ...
        'Menu', ...
        'Modify Parameters','Proceed','Proceed');
    switch choice
        case 'Modify Parameters'
            choice = 1;
        case 'Proceed'
            choice = 0;
    end

    if choice == 1

        prompt = {
            'Steerable_Sigma :'
            };

        definput = {
            num2str(Steerable_Sigma)
            };
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        
        Steerable_Sigma = answer(1,1); 
        
        close
    end

    if choice == 0
        
        close
    end    
    
end



% #########################################################################
%                   SAVE ALL PARAMETERS NEEDED FOR NEXT SCRIPT
% #########################################################################

Parameters = vertcat(...
    GBlur_Sigma,...
    Mask_Thresh,...
    ROI_Size,ROI_Thresh,...
    Tophat_Sigma,...
    Steerable_Sigma);

Parameters = array2table(Parameters,'RowNames',...
    {'GBlur_Sigma',...
    'Mask_Thresh',...
    'ROI_Size','ROI_Thresh',...
    'Tophat_Sigma',...
    'Steerable_Sigma'});


clearvars -except...
    RawPath RootPath RawName...
    Raw Raw_BGSub Raw_Mask Raw_ROIsMask Raw_Tophat... 
    Parameters...
    AngleMap...
    nGridX nGridY...
    pixSize
