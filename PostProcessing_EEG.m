%**************************************************************************
%******************* ERP Features Extraction per channel/per region *******
%**************************************************************************
%*****4-Extract statistics Features, Power Spectral Density of Channels****
%******5- Latency, Amplitude, Area, ***************************************
%*******6- Empirical Mode Decomposition: IMF 1 to IMF 3********************
%*******7- Energy, Entropy, Statistical features of IMFs*******************
%***************************written by Lina Abou-Abbas*********************
%***************************Email: lina.abou-abbas@mail.mcgill.ca**********

clc
% close all
clear all
BaselineCorrect=-100;%in ms
%epoch_limit_start=-0.100;%in s
%epoch_limit_end=+1.700;%in s

Folder_name=[int2str(epoch_limit_start),'ms-',int2str(epoch_limit_end),'ms'];

% Channel_location_file='C:\Users\linaa\OneDrive - McGill University\Washington Data for Analysis\misc\GSN129.sfp';
% Channel_location_file_1020='C:\Users\linaa\OneDrive - McGill University\Washington Data for Analysis\misc\standard_1020_bucanl19.elc';
% Event_Code_file='C:\Users\linaa\OneDrive - McGill University\Washington Data for Analysis\misc\EventCode-mod.txt';
% check_path='C:\Users\linaa\OneDrive - McGill University\Washington Data for Analysis\qcr';
% study_path='C:\Users\linaa\OneDrive - McGill University\Washington Data for Analysis\Study\ERP';

Channel_location_file='C:\Users\linaa\OneDrive - McGill University\London Data for Analysis\misc\GSN-HydroCel-129_EEGLAB.sfp';
Channel_location_file_1020='C:\Users\linaa\OneDrive - McGill University\London Data for Analysis\misc\standard_1020_bucanl19.elc';
Event_Code_file='C:\Users\linaa\OneDrive - McGill University\London Data for Analysis\misc\EventCode-mod.txt';
check_path='C:\Users\linaa\OneDrive - McGill University\London Data for Analysis\qcr';
study_path='C:\Users\linaa\OneDrive - McGill University\London Data for Analysis\Study\ERP';



check_dir=dir(check_path); % list of folders in current path
sr=length(check_dir);%number of folders in current path
nb_sig_av=0;

for sn=3:3%sr
    
    try;
    
        %sr % the 8 first two elements in input path struct do not contain names of folders
        check_dir(sn).name
        mkdir(study_path,check_dir(sn).name);
        %ses-m06 is for the folder that contains eeg data for phase I- of 6 months
        %old. for the folder of 12 months use 'ses-m12' instead of ses-m06
        study_erp_path=fullfile([study_path],check_dir(sn).name,'ses-m06\');
        mkdir(study_erp_path);
        % ERP_path=fullfile([study_set_fdt_path],'ERP')
        set_fdt_path=fullfile([check_path],check_dir(sn).name,'\ses-m06\eeg')
        set_fdt_dir=dir([set_fdt_path])
        %set_fdt_dir=dir([check_path],check_dir(sn).name)
        signal_set_name=set_fdt_dir(4).name;
        c=char(signal_set_name);
        signal_set_name_no_ext=c(1:end-12);
        signal_set_path=fullfile([check_path,check_dir(sn).name],signal_set_name);
        EEG = pop_loadset('filename',signal_set_name,'filepath',set_fdt_path);
        EEG = eeg_checkset( EEG );

           for i=1:length(EEG.event)
              if strcmp(EEG.event(i).type, 'Face') % if strcmp(EEG.event(i).type, 'st+1') 
              cont=1; break
               end
           end

   %% 
   
% Removed flagged channels and time segments

    if cont==1
    nb_sig_av=nb_sig_av + 1    
% % To add this section with Washington data
%             for i=1:length(EEG.event);
%                 if strcmp(EEG.event(i).type, 'badt')
%                 EEG.event(i).type = 'boundary';
%                 end
%             end
            sprintf('%s','Purging flagged channels...\n');
            try
            EEG = pop_marks_select_data(EEG,'channel marks',[],'labels',{'manual'},'remove','on');
            end
            try
            EEG = pop_marks_select_data(EEG,'time marks',[],'labels',{'manual'},'remove','on');
            end
            try
            EEG = pop_marks_select_data(EEG,'component marks',[],'labels',{'manual'},'remove','on');
            end
            EEG = eeg_checkset(EEG);
            % EEG = pop_saveset( EEG, 'filename',[signal_set_name_no_ext,'manual_removed.set'], ...
            %    'filepath',study_erp_path);
            %purge unnecessary fields...
             for i=1:length(EEG.marks.time_info);
                 EEG.marks.time_info(i).flags=[];
             end

            EEG.data_sd_ep=[];
            EEG.c_data_sd_ep=[];
            EEG.data_sd_ch=[];
            EEG.c_data_sd_ch=[];
            EEG.m_neigbr_r_ch=[];
            EEG.c_neigbr_r_ch=[];
            EEG.m_neigbr_r_lat=[];
            EEG.c_neigbr_r_lat=[];
            EEG.amica=[];
            EEG.icaact_sd1_lat=[];
            EEG.c_icaact1_sd_lat=[];
            EEG.icaact_sd2_lat=[];
            EEG.c_icaact_sd2_lat=[];
            EEG.icaact_ta_lat=[];
            EEG.c_icaact_ta_lat=[];
            EEG.icaact_b_lat=[];
            EEG.c_icaact_b_lat=[];
            EEG.icaact = [];
            EEG.icawinv = [];
            EEG.icasphere = [];
            EEG.icaweights = [];
            EEG.icachansind = [];
            % EEG.chaninfo.nosedir='-X';
            EEG=eeg_checkset(EEG);
            
            % EEG = pop_saveset( EEG, 'filename',[signal_set_name_no_ext,'Rmv_bad_Marks.set'], ...
            %   'filepath',[study_path,'\',check_dir(sn).name]);%save EEG after removing bad marks in the study path for comparison purpose
%             manual_chans = find(EEG.marks.chan_info(1).flags)'; % save flagged channel indices into variable
%             EEG = eeg_interp(EEG,manual_chans,'spherical');%interpolate bad channels
            tmpEEG=EEG;
            %interpolate
             EEG = warp_locs(EEG,Channel_location_file,'transform', ...
                            [0,0,0,0,0,-1.57,1,1,1],'manual','off');

            EEG = interp_mont( EEG,Channel_location_file_1020,'manual','off');
            EEG = pop_reref(EEG,[]);
            
            EEG = pop_eegfiltnew(EEG,[], 1,[],1,[],0);
            EEG = pop_eegfiltnew(EEG,[],30,[],0,[],0);

            %% 

            tmpEEG=EEG;
            %% 

            % epoch Face 
%             try;
%                 EEG = pop_epoch( tmpEEG, { 'Face' }, ... %st+1
%                                 [-0.2 0.8], 'newname', 'Face', 'epochinfo', 'yes');
%                 EEG = pop_rmbase( EEG, [-200    0]);
%                 EEG.condition = 'Face';
% %                 EEG = pop_saveset( EEG, 'filename',[signal_set_name_no_ext,'_face_seg.set'], ...
% %                'filepath',study_erp_path);
%             catch;
%                 disp(['SKIPPING face condition for ' EEG.filename ': No segments remaining.']);
%             end;
%             EEG_epoch_face=EEG
%             EEG=tmpEEG;
%             % epoch Object
% %             try;
% %                 EEG = pop_epoch( tmpEEG, { 'st+2' }, ...
% %                                 [-0.2 0.8], 'newname', 'Object', 'epochinfo', 'yes');
% %                 EEG = pop_rmbase( EEG, [-200    0]);
% %                 EEG.condition = 'Object';
% %                 EEG = pop_saveset( EEG, 'filename',[signal_set_name_no_ext,'_object_seg.set'], ...
% %                'filepath',study_erp_path);
% %             catch;
% %                 disp(['SKIPPING Object condition for ' EEG.filename ': No segments remaining.']);
% %             end;
% %             EEG_epoch_object=EEG
%             EEG=tmpEEG;
%             
%             % epoch Noise
%             try;
%                 EEG = pop_epoch( tmpEEG, { 'Nois' }, ...
%                                 [-0.2 0.8], 'newname', 'Noise', 'epochinfo', 'yes');
%                 EEG = pop_rmbase( EEG, [-200    0]);
%                 EEG.condition = 'Noise';
% %                 EEG = pop_saveset( EEG, 'filename',[signal_set_name_no_ext,'_noise_seg.set'], ...
% %                'filepath',study_erp_path);
%             catch;
%                 disp(['SKIPPING noise condition for ' EEG.filename ': No segments remaining.']);
%             end;
% 
%             EEG_epoch_noise=EEG
%             EEG=tmpEEG;

            % EEG=pop_reref(EEG,[]);
            % Creating the event list using a defined file called EyeGazeEventCode, using a fully specified path

            % A text file with the event list saved in elist.txt in the same path
            EEG = pop_editeventlist(EEG, 'List',Event_Code_file,'ExportEL','EventList.txt', ...
                'BoundaryNumeric', {-99});


            % Copy event labels into the EEG structure
            EEG = pop_overwritevent( EEG, 'codelabel');

            % Create Bin-Based EEG Epochs
            % 'pre' means use the prestimulus period for baseline correction
            EEG = pop_epochbin( EEG , [-200.0  800.0], [-200 0])
            EEG=pop_rmbase(EEG,[-200 0]);
            % Saving after Epoching
            % Setname is S1_EEG_elist_be
            %EEG.setname=[signal_set_name_no_ext,'_NewEventCode']
            %EEG=pop_saveset(EEG, 'filename', [EEG.setname,'.set'],'filepath',set_fdt_path)%'[study_path,check_dir(sn).name]);
            % Compute Averaged ERP: classic(time domain), Spectral Domain: Total Power Spectrum ('TFFT'),the Evoked Power Spectrum ('EFFT')
            ERP = pop_averager( EEG , 'Criterion', 'all')%, 'SEM', 'on','Compute','ERP','ExcludeBoundary','on');
            % ERP_TFFT = pop_averager( EEG , 'Compute', 'TFFT', 'Criterion', 'all', 'ExcludeBoundary', 'on', 'SEM', 'on','TaperWindow');
            % ERP_EFFT = pop_averager( EEG , 'Compute', 'EFFT', 'Criterion', 'all', 'ExcludeBoundary', 'on', 'SEM', 'on','TaperWindow');
            % Setname  is followed by _ERPs
            % Filename is followed by _ERPs.erp
            %% Filter ERP
                        % Channels; No high-pass;
                        % Lowpass cutoff at 30 Hz; Order of the filter = 2.
                        % Type of filter = "Butterworth"; Do not remove DC offset
%             ERP = pop_filterp( ERP,1:ERP.nchan, 'Cutoff',80, 'Design', 'butter', 'Filter', 'lowpass', 'Order',2 );

            % Note that you will need to replace the path with the actual location in your file system
            % ERP=pop_savemyerp(ERP, 'erpname', 'test_ERP','filename', 'test_ERP.erp', 'filepath', study_erp_path);
%             ERP=pop_savemyerp(ERP, 'erpname', [signal_set_name_no_ext,'_ERP'],'filename', ...
%                 [signal_set_name_no_ext,'_ERP.erp'] , 'filepath', study_erp_path);

           

%             save([study_erp_path,'\',signal_set_name_no_ext,'_ERP_face.mat'],'ERP')
%             save([study_erp_path,'\',signal_set_name_no_ext,'_ERP_object.mat'],'ERP_object')





            %% 
            %This section is to compute the peaks of ERP signal, to extract the value
            %and the position of the maximum positive peak and the minimum negative
            %peak

            % ERP per channel
            [ERP,ERP_area] = pop_geterpvalues( ERP, [ 0 789],  ERP.nbin ,  1:ERP.nchan , 'Baseline', 'pre', 'Binlabel',...
                'on', 'FileFormat', 'long', 'Filename', 'text.txt', 'Fracreplace', 'NaN', 'InterpFactor',1,...
                'Measure', 'areat', 'PeakOnset',1, 'Resolution',3, 'SendtoWorkspace', 'off' );
            
            
            [ERP,ERP_Pk_amp_pos,ERP_Pk_lat_pos] = pop_geterpvalues( ERP, [ 0 789],  ERP.nbin ,  1:ERP.nchan ,...
                'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
             'IncludeLat', 'yes', 'InterpFactor',  1,  'Measure', 'peaklatbl','Measure', 'peakampbl',...
             'Neighborhood',  3, 'PeakOnset',  1, 'Peakpolarity', 'positive', 'Peakreplace',...
             'absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
            ERP_Pk_lat_pos=cell2mat(ERP_Pk_lat_pos{1,:});
            
             [ERP,ERP_Pk_amp_neg,ERP_Pk_lat_neg] = pop_geterpvalues( ERP, [ 0 789],  ERP.nbin ,  1:ERP.nchan ,...
                'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
             'IncludeLat', 'yes', 'InterpFactor',  1, 'Measure', 'peaklatbl', 'Measure', 'peakampbl',...
             'Neighborhood',  3, 'PeakOnset',  1, 'Peakpolarity', 'negative', 'Peakreplace',...
             'absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
            ERP_Pk_lat_neg=cell2mat(ERP_Pk_lat_neg{1,:});
            % 
            % %This section is to compute P170, N290 and P400 of the ERP signals
            % %P400=ERP_MEASURES;
            % %delta=1-4, theta=4-8, alpha=8-13, beta=13-30, gamma=30-80
            % %position of peaks based on 
            % %% P1 for gaze shift intervals different than other bins
            [ERP,ERP_Pk_amp_face_P1,ERP_Pk_lat_face_P1] = pop_geterpvalues( ERP, [ 100 219],  ERP.nbin,  1:ERP.nchan ,...
                'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
             'IncludeLat', 'yes', 'InterpFactor',  1, 'Measure', 'peaklatbl','Measure', 'peakampbl', ...
             'Neighborhood',  3, 'PeakOnset',  1, 'Peakreplace','absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
            ERP_Pk_lat_face_P1=cell2mat(ERP_Pk_lat_face_P1{1,:});
            % 
            % 
            % % Earliest N1 for gaze shift
            % 
            [ERP,ERP_Pk_amp_face_N290,ERP_Pk_lat_face_N290] = pop_geterpvalues( ERP, [ 220 319], ERP.nbin,  1:ERP.nchan ,...
                'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
             'IncludeLat', 'yes', 'InterpFactor',  1, 'Measure', 'peaklatbl','Measure', 'peakampbl',...
             'Neighborhood',  3, 'PeakOnset',  1, 'Peakreplace','absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
            ERP_Pk_lat_face_N290=cell2mat(ERP_Pk_lat_face_N290{1,:});
            % 
             [ERP,ERP_Pk_amp_face_P400,ERP_Pk_lat_face_P400] = pop_geterpvalues( ERP, [ 320 540],  ERP.nbin,  1:ERP.nchan ,...
                  'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN', 'IncludeLat', ...
                  'yes', 'InterpFactor',  1, 'Measure', 'peakampbl', 'Neighborhood',  3, 'PeakOnset',  1,  ...
                  'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
             ERP_Pk_lat_face_P400=cell2mat(ERP_Pk_lat_face_P400{1,:});
              
             [ERP,ERP_Pk_amp_face_LPC1,ERP_Pk_lat_face_LPC1] = pop_geterpvalues( ERP, [ 541 789],  ERP.nbin,  1:ERP.nchan , ...
                  'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
                  'IncludeLat', 'yes', 'InterpFactor',  1, 'Measure', 'peaklatbl', 'Measure', 'peakampbl',...
                  'Neighborhood',  3, 'PeakOnset',  1, 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
             ERP_Pk_lat_face_LPC1=cell2mat(ERP_Pk_lat_face_LPC1{1,:});
             
             ERP_features_chan_onerow=[ERP_area, ERP_Pk_amp_pos,ERP_Pk_lat_pos,ERP_Pk_amp_neg,ERP_Pk_lat_neg,...
                 ERP_Pk_amp_face_P1,ERP_Pk_lat_face_P1,ERP_Pk_amp_face_N290,ERP_Pk_lat_face_N290,...
                 ERP_Pk_amp_face_P400,ERP_Pk_lat_face_P400,ERP_Pk_amp_face_LPC1,ERP_Pk_lat_face_LPC1];
             
           save([study_erp_path,'\',signal_set_name_no_ext,'_Allfeatures_ERP_onevector.mat'],'ERP_features_chan_onerow');

           %% 

             %London epochs times are from -200 to 800ms Washington epochs from
%-100 to 1600ms -----> to merge databases in face stimulus we consider only
% ERP components from 0 to 800 ms
              %            [ERP,ERP_Pk_amp_face_LPC2,ERP_Pk_lat_face_LPC2] = pop_geterpvalues( ERP, [ 790 900],  ERP.nbin,  1:ERP.nchan , ...
%                   'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
%                   'IncludeLat', 'yes', 'InterpFactor',  1, 'Measure', 'peaklatbl', 'Measure', 'peakampbl',...
%                   'Neighborhood',  3, 'PeakOnset',  1, 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
%               
            for i=1:EEG.nbchans
                for j=1:EEG.pnts
                    avg_EEG_face(i,j)=mean(EEG.data(i,j,:));
                end
            end
            % 
            for i=1:EEG.nbchans
                
            [PSD_dB_avg(i,:),freqs] =spectopo(ERP.bindata(i,:),0,EEG.srate,'freqrange',[1 80],...
                'plot','off','logtrials', 'off','overlap',250,'wintype','hamming');  
                           
            end
            
            
            deltaIndx = find(freqs>1 & freqs<4);
            thetaIndx = find(freqs>4 & freqs<8);
            alphaIndx = find(freqs>8 & freqs<13);
            betaIndx  = find(freqs>13 & freqs<30);
            gammaIndx = find(freqs>30 & freqs<50);
            
              for i=1:EEG.nbchans
                
            %Compute absolute power per each channel per each frequency band 1st method
            AbsdeltaPower(i) = mean(10.^(PSD_dB_avg(i,deltaIndx)/10));
            AbsthetaPower(i) = mean(10.^(PSD_dB_avg(i,thetaIndx)/10));
            AbsalphaPower(i) = mean(10.^(PSD_dB_avg(i,alphaIndx)/10));
            AbsbetaPower (i) = mean(10.^(PSD_dB_avg(i,betaIndx)/10));
            AbsgammaPower(i) = mean(10.^(PSD_dB_avg(i,gammaIndx)/10));
           
            TotalPower(i)= AbsdeltaPower(i)+AbsthetaPower(i)+ AbsalphaPower(i) +AbsbetaPower(i)+AbsgammaPower(i);

            ReldeltaPower(i)=AbsdeltaPower(i)/TotalPower(i);
            RelthetaPower(i)=AbsthetaPower(i)/TotalPower(i);
            RelalphaPower(i)=AbsalphaPower(i)/TotalPower(i);
            RelbetaPower(i)=AbsbetaPower(i)/TotalPower(i);
            RelgammaPower(i)=AbsgammaPower(i)/TotalPower(i);
            
              end
              
            Power_features_chan_onerow=[ AbsdeltaPower, AbsthetaPower, AbsalphaPower ,AbsbetaPower,AbsgammaPower, ...
                ReldeltaPower, RelthetaPower, RelalphaPower ,RelbetaPower,RelgammaPower];           
    
           save([study_erp_path,'\',signal_set_name_no_ext,'_Allfeatures_Power_onevector.mat'],'Power_features_chan_onerow');

            % 
            %% %% ******************** Extract IMF and Calculate Statistical Features ********************
            %% %% Total of 18 features per bin per channel
            for chan=1:ERP.nchan
%                IMF_features_chan=[]
%                 for m=1:ERP.nbin
                   IMF_features_vec=[]

                   imf_ERP=emd(ERP.bindata(chan,:),'SiftRelativeTolerance',0.0001); 
                        for i=1:3 % first three imf_num
                            % ENERGY
                            imf_ener=sum(imf_ERP(:,i).^2);
                            % ENTROPY (Shannon)
                            temp = double(imf_ERP(:,i));
                            temp = temp(temp>0).^2;
                            imf_entr = -sum(temp.*log(eps+temp));
                            % MEAN
                            imf_mean=mean(imf_ERP(:,i));
                            % STANDARD DEVIATION
                            imf_std=std(imf_ERP(:,i));
                            % THIRD MOMENT
                            imf_mom=moment(imf_ERP(:,i),3);
                            %         disp(int2str(i))
                            imf_skew=skewness(imf_ERP(:,i));
                            %         disp(int2str(i))
                            IMF_features_vec=[IMF_features_vec,imf_ener,imf_entr,imf_mean,imf_std,imf_mom,imf_skew];
                        end

              IMF_features_chan(chan,:)=IMF_features_vec;
            end
            IMF_features_chan_onerow=[];
            [i,j]=size(IMF_features_chan)
            n=1;
            for k=1:i 
                for l=1:j 
                   IMF_features_chan_onerow(1,n)= IMF_features_chan(k,l);
                   n=n+1;
                end
            end

           save([study_erp_path,'\',signal_set_name_no_ext,'_Allfeatures_IMF_onevector.mat'],'IMF_features_chan_onerow');
          
           
           %%  %%    % ERP_sys1020_region
           
                       %% %% Averaging ERP to calculate ERP per region per bin if 10/20 system is used
            %Frontal:AF8,Fpz,AF7,F3,Fz,F4:1,2,3,4,5,6
            %Central:C3,Cz,C4:9,10,11
            %Occipital, Parietal:P3,Pz,P4,PO7,PO8, Oz: 14,15,16,17,18,19
            %Temporal:FT8,TP8,FT7,TP7: 7,8,12,13

            ERP_sys1020_region = pop_erpchanoperator( ERP, {  'nch1 = (ch1 + ch2 + ch3 + ch4 + ch5+ ch6)/6 label Front',...
              'nch2=( ch9+ ch10+ ch11)/3 label Cent',...
              'nch3= (ch7 + ch8 + ch12 + ch13)/4 label Temp',...
              'nch4= (ch14 + ch15 + ch16+ ch17 + ch18+ ch19)/6 label Occ_Par'} ,...
             'ErrorMsg', 'popup', 'KeepLocations',  0, 'Warning', 'on' );

 
%             ERP_sys1020_region=pop_savemyerp(ERP_sys1020_region, 'erpname', [signal_set_name_no_ext,'_ERP_sys1020_4regions.set'],'filename', ...
%                 [signal_set_name_no_ext,'_ERP_sys1020_4regions.erp'] , 'filepath', study_erp_path);
% 
%                 
                    
            [ERP_sys1020_region,ERP_sys1020_region_area] = pop_geterpvalues( ERP_sys1020_region, [ 0 789],  ERP_sys1020_region.nbin ,  1:ERP_sys1020_region.nchan , 'Baseline', 'pre', 'Binlabel',...
                'on', 'FileFormat', 'long', 'Filename', 'text.txt', 'Fracreplace', 'NaN', 'InterpFactor',1,...
                'Measure', 'areat', 'PeakOnset',1, 'Resolution',3, 'SendtoWorkspace', 'off' );
            
            
            [ERP_sys1020_region,ERP_sys1020_region_Pk_amp_pos,ERP_sys1020_region_Pk_lat_pos] = pop_geterpvalues( ERP_sys1020_region, [ 0 789],  ERP_sys1020_region.nbin ,  1:ERP_sys1020_region.nchan ,...
                'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
             'IncludeLat', 'yes', 'InterpFactor',  1,  'Measure', 'peaklatbl','Measure', 'peakampbl',...
             'Neighborhood',  3, 'PeakOnset',  1, 'Peakpolarity', 'positive', 'Peakreplace',...
             'absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
            ERP_sys1020_region_Pk_lat_pos=cell2mat(ERP_sys1020_region_Pk_lat_pos{1,:});
            
             [ERP_sys1020_region,ERP_sys1020_region_Pk_amp_neg,ERP_sys1020_region_Pk_lat_neg] = pop_geterpvalues( ERP_sys1020_region, [ 0 789],  ERP_sys1020_region.nbin ,  1:ERP_sys1020_region.nchan ,...
                'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
             'IncludeLat', 'yes', 'InterpFactor',  1, 'Measure', 'peaklatbl', 'Measure', 'peakampbl',...
             'Neighborhood',  3, 'PeakOnset',  1, 'Peakpolarity', 'negative', 'Peakreplace',...
             'absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
            ERP_sys1020_region_Pk_lat_neg=cell2mat(ERP_sys1020_region_Pk_lat_neg{1,:});
          
            [ERP_sys1020_region,ERP_sys1020_region_Pk_amp_face_P1,ERP_sys1020_region_Pk_lat_face_P1] = pop_geterpvalues( ERP_sys1020_region, [ 100 219],  ERP_sys1020_region.nbin,  1:ERP_sys1020_region.nchan ,...
                'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
             'IncludeLat', 'yes', 'InterpFactor',  1, 'Measure', 'peaklatbl','Measure', 'peakampbl', ...
             'Neighborhood',  3, 'PeakOnset',  1, 'Peakreplace','absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
            ERP_sys1020_region_Pk_lat_face_P1=cell2mat(ERP_sys1020_region_Pk_lat_face_P1{1,:});
            % 
            % 
            % % Earliest N1 for gaze shift
            % 
            [ERP_sys1020_region,ERP_sys1020_region_Pk_amp_face_N290,ERP_sys1020_region_Pk_lat_face_N290] = pop_geterpvalues( ERP_sys1020_region, [ 220 319], ERP_sys1020_region.nbin,  1:ERP_sys1020_region.nchan ,...
                'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
             'IncludeLat', 'yes', 'InterpFactor',  1, 'Measure', 'peaklatbl','Measure', 'peakampbl',...
             'Neighborhood',  3, 'PeakOnset',  1, 'Peakreplace','absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
            ERP_sys1020_region_Pk_lat_face_N290=cell2mat(ERP_sys1020_region_Pk_lat_face_N290{1,:});
            % 
             [ERP_sys1020_region,ERP_sys1020_region_Pk_amp_face_P400,ERP_sys1020_region_Pk_lat_face_P400] = pop_geterpvalues( ERP_sys1020_region, [ 320 540],  ERP_sys1020_region.nbin,  1:ERP_sys1020_region.nchan ,...
                  'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN', 'IncludeLat', ...
                  'yes', 'InterpFactor',  1, 'Measure', 'peakampbl', 'Neighborhood',  3, 'PeakOnset',  1,  ...
                  'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
             ERP_sys1020_region_Pk_lat_face_P400=cell2mat(ERP_sys1020_region_Pk_lat_face_P400{1,:});
              
             [ERP_sys1020_region,ERP_sys1020_region_Pk_amp_face_LPC1,ERP_sys1020_region_Pk_lat_face_LPC1] = pop_geterpvalues( ERP_sys1020_region, [ 541 789],  ERP_sys1020_region.nbin,  1:ERP_sys1020_region.nchan , ...
                  'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'long', 'Fracreplace', 'NaN',...
                  'IncludeLat', 'yes', 'InterpFactor',  1, 'Measure', 'peaklatbl', 'Measure', 'peakampbl',...
                  'Neighborhood',  3, 'PeakOnset',  1, 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'off' );
             ERP_sys1020_region_Pk_lat_face_LPC1=cell2mat(ERP_sys1020_region_Pk_lat_face_LPC1{1,:});
             
             ERP_sys1020_region_features_chan_onerow=[ERP_sys1020_region_area, ERP_sys1020_region_Pk_amp_pos,ERP_sys1020_region_Pk_lat_pos,ERP_sys1020_region_Pk_amp_neg,ERP_sys1020_region_Pk_lat_neg,...
                 ERP_sys1020_region_Pk_amp_face_P1,ERP_sys1020_region_Pk_lat_face_P1,ERP_sys1020_region_Pk_amp_face_N290,ERP_sys1020_region_Pk_lat_face_N290,...
                 ERP_sys1020_region_Pk_amp_face_P400,ERP_sys1020_region_Pk_lat_face_P400,ERP_sys1020_region_Pk_amp_face_LPC1,ERP_sys1020_region_Pk_lat_face_LPC1];
             
             
           save([study_erp_path,'\',signal_set_name_no_ext,'_Allfeatures_ERP_region_onevector.mat'],'ERP_sys1020_region_features_chan_onerow');
           %% 

             for i=1:ERP_sys1020_region.nchan
                
            [PSD_dB_region(i,:),freqs] =spectopo(ERP_sys1020_region.bindata(i,:),0,EEG.srate,'freqrange',[1 80],...
                'plot','off','logtrials', 'off','overlap',250,'wintype','hamming');  
             end              
            
            
            
            deltaIndx = find(freqs>1 & freqs<4);
            thetaIndx = find(freqs>4 & freqs<8);
            alphaIndx = find(freqs>8 & freqs<13);
            betaIndx  = find(freqs>13 & freqs<30);
            gammaIndx = find(freqs>30 & freqs<50);
            
              for i=1:ERP_sys1020_region.nchan
                
            %Compute absolute power per each channel per each frequency band 1st method
            AbsdeltaPower_region(i) = mean(10.^(PSD_dB_region(i,deltaIndx)/10));
            AbsthetaPower_region(i) = mean(10.^(PSD_dB_region(i,thetaIndx)/10));
            AbsalphaPower_region(i) = mean(10.^(PSD_dB_region(i,alphaIndx)/10));
            AbsbetaPower_region(i) = mean(10.^(PSD_dB_region(i,betaIndx)/10));
            AbsgammaPower_region(i) = mean(10.^(PSD_dB_region(i,gammaIndx)/10));
           
            TotalPower_region(i)= AbsdeltaPower_region(i)+AbsthetaPower_region(i)+ AbsalphaPower_region(i) +AbsbetaPower_region(i)+AbsgammaPower_region(i);

            ReldeltaPower_region(i)=AbsdeltaPower_region(i)/TotalPower_region(i);
            RelthetaPower_region(i)=AbsthetaPower_region(i)/TotalPower_region(i);
            RelalphaPower_region(i)=AbsalphaPower_region(i)/TotalPower_region(i);
            RelbetaPower_region(i)=AbsbetaPower_region(i)/TotalPower_region(i);
            RelgammaPower_region(i)=AbsgammaPower_region(i)/TotalPower_region(i);
            
              end
              
            Power_region_features_chan_onerow=[ AbsdeltaPower_region, AbsthetaPower_region, AbsalphaPower_region ,AbsbetaPower_region,AbsgammaPower_region, ...
                ReldeltaPower_region, RelthetaPower_region, RelalphaPower_region ,RelbetaPower_region,RelgammaPower_region];           
    
           save([study_erp_path,'\',signal_set_name_no_ext,'_Allfeatures_Power_region_onevector.mat'],'Power_region_features_chan_onerow');

            %% %% ******************** Extract IMF per region and Calculate Statistical Features ********************
            %% %% Total of 18 features per bin per channel
            for chan=1:ERP_sys1020_region.nchan
%                IMF_features_chan=[]
%                 for m=1:ERP.nbin
                   IMF_features_region_vec=[]

                   imf_ERP_region=emd(ERP_sys1020_region.bindata(chan,:),'SiftRelativeTolerance',0.0001); 
                        for i=1:3 % first three imf_num
                            % ENERGY
                            imf_ener_region=sum(imf_ERP_region(:,i).^2);
                            % ENTROPY (Shannon)
                            temp = double(imf_ERP_region(:,i));
                            temp = temp(temp>0).^2;
                            imf_entr_region = -sum(temp.*log(eps+temp));
                            % MEAN
                            imf_mean_region=mean(imf_ERP_region(:,i));
                            % STANDARD DEVIATION
                            imf_std_region=std(imf_ERP_region(:,i));
                            % THIRD MOMENT
                            imf_mom_region=moment(imf_ERP_region(:,i),3);
                            %         disp(int2str(i))
                            imf_skew_region=skewness(imf_ERP_region(:,i));
                            %         disp(int2str(i))
                            IMF_features_region_vec=[IMF_features_region_vec,imf_ener_region,imf_entr_region,imf_mean_region,imf_std_region,imf_mom_region,imf_skew_region];
                        end

              IMF_features_region(chan,:)=IMF_features_region_vec;
            end
            IMF_features_region_onerow=[];
            [i,j]=size(IMF_features_region)
            n=1;
            for k=1:i 
                for l=1:j 
                   IMF_features_region_onerow(1,n)= IMF_features_region(k,l);
                   n=n+1;
                end
            end

           save([study_erp_path,'\',signal_set_name_no_ext,'_Allfeatures_IMF_region_onevector.mat'],'IMF_features_region_onerow');
           nb_good_trials=EEG.trials
           csvwrite([study_erp_path,'\',signal_set_name_no_ext,'_nb_good_trials.csv'],nb_good_trials);
           
    end
    
    
    catch; 
    
    end;
sn
end
