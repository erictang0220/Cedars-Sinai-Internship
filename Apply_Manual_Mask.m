%% Creat Mask
% Creat mask if original mask doesn't exist

MaskView = 'Axis'; % This is the orientation of the matrix. you can start with
FlagMask = true; % meaning use mask
IsManual_Mask = true; % meaning draw manual mask.
%RawdataFT = UNIC_B0Map.complexmap;

n=1;
AllPhasemap(1).Mask=ones(size(CTmag)); % UNIC_B0Map.mag
% AllPhasemap(1).Mask=UNIC_B0Map.Parameters.Mask; % to generate a dummy mask for the code to start with.
AllPhasemap(n).ManualMask=0; % 

if ~exist('Manual_Mask.mat')
    mag_all = CTmag.*AllPhasemap(1).Mask;% UNIC_B0Map.mag.*AllPhasemap(1).Mask; %sqrt(sum(abs(permute(RawdataFT(:,:,:,:,1),[1 2 3 5 4])),5)).*AllPhasemap(1).Mask;
    mag_all=mag_all/(max(mag_all(:)));
    % Readout FOV at center
    RD_FOV=[(round(size(mag_all,1)/2)+1-round(size(mag_all,1)/4))...
        :(round(size(mag_all,1)/2)+1+round(size(mag_all,1)/4))];
    % Change the view angle of rough mask
    Rawmask=AllPhasemap(1).Mask;
    clear Newmask Mask_m; %what does this do
    %%
    switch MaskView
        case 'Axis'
            centMask=zeros(size(mag_all));
            centMask(RD_FOV,1:end,1:end,1:end)=1;
            Rawmask=AllPhasemap(n).Mask;
        case 'Coronal'
            mag_all=permute(mag_all, [3 1 2]);% get the side view and readout to be dim 2
            centMask=zeros(size(mag_all));
            centMask(1:end, RD_FOV,1:end,1:end)=1;
            Rawmask=permute(AllPhasemap(n).Mask,[3 1 2]);
        case 'Sagittal'
            mag_all=permute(mag_all, [2 3 1]);% get the side view and readout to be dim 2
            centMask=zeros(size(mag_all));
            centMask(1:end,1:end,RD_FOV,1:end)=1;
            Rawmask=permute(AllPhasemap(n).Mask,[2 3 1]);
        otherwise
            centMask=zeros(size(mag_all));
            centMask(RD_FOV,1:end,1:end,1:end)=1;
            Rawmask=AllPhasemap(n).Mask;
    end
    %%
    mag_all=mag_all.*centMask;
    Generatecontour=1;
    % create mask
    Mask_m=nan;
    for ms=1:2:(size(mag_all,3)-1);
        mag=mag_all(:,:,ms);
        if FlagMask
            if (IsManual_Mask && AllPhasemap(n).ManualMask==0)
                if n==1;%only active when analyzing the first file
                    %Charles 27Jun2019, add to show slice# as well
                    disp(['Mask ', num2str(ms/size(mag_all,3)*100),'%, ', num2str(ms), '/', num2str(size(mag_all,3)), ' slice'])% to display in a better orientation
                    nmag=mag/max(mag(:))-median(mag(:));%for imshow normalize the figure
                    
                  if sum(sum(Rawmask(:,:,ms)))>0
                    if isnan(Mask_m)
                        imagesc(nmag);colormap gray;
                        set(gcf, 'units', 'normalized', 'Position', [0.045, 0.3, 0.4, 0.5]) %Charles 30Jan2019, show the figure at a better location
                        %Orginally: set(gcf, 'Position', [100, 100, 1000, 1000])
                        [Mask_m tmp_xi tmp_yi]=roipoly;
                        0
                    else
                        imagesc(nmag);colormap gray;
                        set(gcf, 'units', 'normalized', 'Position', [0.045, 0.3, 0.4, 0.5]) %Charles 30Jan2019, show the figure at a better location
                        %Orginally: set(gcf, 'Position', [100, 100, 1000, 1000])
                        h = impoly(gca,[tmp_xi tmp_yi]); 
                        Generatecontour=input('Re-contour?');
                        Mask_m = createMask(h); %BW contains the mask which you just altered.
                        pos = getPosition(h); 
                        tmp_xi=pos(:,1);
                        tmp_yi=pos(:,2);
                        if Generatecontour % overwrite mask if need more points
                            imagesc(nmag);colormap gray;
                            [Mask_m tmp_xi tmp_yi]=roipoly;
                        end    
                    end
                  end
                    Mask_m(isnan(Mask_m))=0;
                    Newmask(:,:,ms)=imfill(Rawmask(:,:,ms).*Mask_m,'holes');

                    Newmask(:,:,ms+1)=Newmask(:,:,ms);
                    if ms==(size(Rawmask,3)-1);
                        AllPhasemap(n).ManualMask=1;
                    end
                else
                    AllPhasemap(n).Mask(:,:,ms)= AllPhasemap(1).Mask(:,:,ms);
                    if ms==(size(AllPhasemap(n).Mask,3)-1);
                        AllPhasemap(n).ManualMask=1;
                    end
                end
                
            end
            % mask = repmat(AllPhasemap(n).Mask(:,:,ms),[1,1,size(phase,3)]);
            % phase(~mask)=nan;
            
        end
    end
    

    %fill holes
%     for s=1:size(AllPhasemap(n).Mask,3)
%         Newmask(:,:,s)= imfill(Newmask(:,:,s),'holes');
%     end
    switch MaskView
        case 'Axis'
             AllPhasemap(n).Mask=Newmask;
        case 'Coronal'
            AllPhasemap(n).Mask=permute(Newmask,[2 3 1]); 
        case 'Sagittal'
            AllPhasemap(n).Mask=permute(Newmask,[3 1 2]); 
        otherwise
            centMask=zeros(size(mag_all));
            centMask(RD_FOV,1:end,1:end,1:end)=1;
            Rawmask=AllPhasemap(n).Mask;
    end
    

    Mask=AllPhasemap(n).Mask;
    save('Manual_Mask.mat', 'Mask');
else
    load('Manual_Mask.mat');
    AllPhasemap(n).Mask=Mask;
end
    %clear Mask