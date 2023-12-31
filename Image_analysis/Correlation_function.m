%%  Correlation functions 
% < This script computes correlation function from the linescans using the nematic director >
% Copyright (C) <2023>  <Francine Kolley>

% Francine Kolley 
% Physics of Life, Benjamin M. Friedrich group
% TU_dresden 
% contact: francine.kolley@tu-dresden.de


% #######################################################
% With this script we measure the ACF/CCF for the different combinations
% This script is used after measuring the local nematic order, by using the script nematic_order.m 
% choose the needed channels for the correlation
% #######################################################
% Note that our aim is to compute nomralized and unbiased correlation
% functions and we account for the mean value at the end .
% Unbiased corr. functions lead to a shrinked interval you can believe
% for this make sure that you choose the corr-length long enough to believe
% the peaks you want to measure.
% NO PEAK FILTERING APPLIES HERE!!! WE COMPUTE ALL POSSIBLE CORRELATION
% FUNCTIONS FROM LINESCANS AND CALCULATE THE MEAN OVER THE WHOLE IMAGE!!!


close all

%% Choose the channels 
nChannel1=1; % channel for first substrate
nChannel2=4; %channel for second substrate

%% Meta-informations 
RawPath = strcat(RawPath);
% Open Data
info = imfinfo(RawPath);

%% Read out neccessary parameter for calculation
% get the parameters from Mask
GBlur_Sigma=Parameters{1,1};
ROI_Size=Parameters{3,1};
Tophat_Sigma=Parameters{5,1};

%% Read in intensities
% use same parameters as for the masks
Raw1 = imread(RawPath,nChannel1,'Info',info); %substrate 1
Raw1_BGSub = imgaussfilt(Raw1,GBlur_Sigma); % Gaussian blur #1
Raw1_Tophat=imtophat(Raw1_BGSub,strel('disk',Tophat_Sigma)); %Raw before

Raw2 = imread(RawPath,nChannel2,'Info',info); %substrate 2 
Raw2_BGSub = imgaussfilt(Raw2,GBlur_Sigma); % Gaussian blur #2
Raw2_Tophat=imtophat(Raw2_BGSub,strel('disk',Tophat_Sigma)); %Raw before

%% correlation length and other pre-defintions of values
% Correlation length
corr_length=4; % [\mu m]
crop_length=round(corr_length/(2*pixSize)); %[pixel]

% Rotate Corr2D according to "Angle" 
max_dist=round(crop_length)+10; % +10pixel to corr_length to smaller deviations 
maxlag = 2*max_dist; % Correlation length will be doubled, because of mirror symmetry

ccf_all = nan(nGridY*nGridX,2*maxlag+1); % allocate memory for ACF/CCF

MergedData = cell(nGridY,nGridX);
z=1; % counting number 


%% Computation over the whole grid 
figure(1), hold on
for i=1:nGridY
    for j=1:nGridX
        
        if Raw_ROIsMask(i,j) == 1 % If they are of interest
   
            Angle = AngleMap(i,j); % find right angle for current window
            ccf = zeros(2*maxlag+1,1); % allocate memory for ACF/CCF
            var1 = 0; % reset var1
            var2 = 0; % reset var2
            
            kmax = 10;  % in [pixel] 

            % Do linescans in both, x and y direction and compute from that
            % the 1D correlation function 
            for k=-kmax:2:kmax
                xi = ROI_Size*j+max_dist*[-1 1]*cos((90-Angle)*-1*pi/180)+k*cos((Angle)*pi/180); %x-position
                yi = ROI_Size*i+max_dist*[-1 1]*sin((90-Angle)*-1*pi/180)+k*sin((Angle)*pi/180); %y-position
                linescan1 = improfile( Raw1_Tophat, xi, yi, 2*max_dist+1); %get the profile for substrate 1                           
                linescan1 = double( linescan1 ); 
                linescan2 = improfile( Raw2_Tophat, xi, yi, 2*max_dist+1); %get the profile for substrate 2                                 
                linescan2 = double( linescan2 ); 

                % substract the mean
                linescan1_mean = mean( linescan1 ); 
                linescan1 = linescan1 - linescan1_mean;
                linescan2_mean = mean( linescan2 ); 
                linescan2 = linescan2 - linescan2_mean;
                
                %compute unbiased CCF 
                ccf = ccf + xcorr(linescan1,linescan2,'unbiased')/(2*kmax+1); 

                % variances of individual linescans (Normalization)
                var1 = var1 + var( linescan1 ) / (2*kmax+1); 
                var2 = var2 + var( linescan2 ) / (2*kmax+1);                 
            end 
                                   
            % Normalize CCF by respective standard deviations
            ccf = ccf / sqrt( var1 * var2 ); 
            
            % keep a record
            ccf_all((i-1)*nGridX+j ,:) = ccf;
        
%% plot individual ACFs/CCFs
    maxlag_plot = size( linescan1, 1)-1; 
    lags = (-maxlag:maxlag)*pixSize; % [um] 
    ind = maxlag_plot+1: 2*maxlag_plot+1; % start with 'lag = 0'
    plot( lags(ind), ccf(ind),'Color','#cccaca')
    xlabel('\Delta x [um]')
    ylabel('CCF')
    title('CCF (unbiased, normalized)')

        end % ?                       
    end % j 
end % i

%% Find the valid CCFs/ACFs and calculate the mean
ccf_all_valid= ccf_all(all(~isnan(ccf_all),2),:);  %delete rows with nans
mean_ccf=mean(ccf_all_valid); 
% standart error of the mean
div_factor=size(ccf_all_valid,1);
std_mean_ccf=std(ccf_all_valid,1,'omitnan'); %/sqrt(div_factor);

 plot( lags(ind), mean_ccf(ind),'-','Color','#d13111','LineWidth',1.8) 
 plot(lags(ind),mean_ccf(ind)-std_mean_ccf(ind),'--','Color','#d13111','LineWidth',1.8)
 plot(lags(ind),mean_ccf(ind)+std_mean_ccf(ind),'--','Color','#d13111','LineWidth',1.8)
 ylim([-0.5 1])
 xlim([0 corr_length])
 set(gca,'FontSize',24)



