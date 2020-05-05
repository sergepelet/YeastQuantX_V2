function VarCell = FlatnessMeasure(VarCell)
%Analyses flatness of images based on a buch of images.
% removes high intensity pixels in theses images extract a flatness
% measurement


if isfield(VarCell{1}.Analysis, 'Flatify')
    %Set number of Position averaged to get the flatten image
    NumPos = length(VarCell);
    MaxImage = 50;
    WindowSize = 21;
    if NumPos >MaxImage
        NumPos = MaxImage;
    end
    
    %Get index of which illuminations have to be flatten
    ToFlatten = find(strcmpi(VarCell{1}.Analysis.Flatify, 'yes')== 1);
    if ~isempty(ToFlatten)
        %Loop through all illuminations to flatten
        for I = ToFlatten
            %Get Illum name
            Illum = VarCell{1}.Analysis.LoadIllum{I};
            fprintf(['Extract flatness from ', Illum,'\n'])
            %Get filter Number to access file name
            FilterNum = find(strcmp(VarCell{1}.Experiment.Filter, Illum)==1);
            
            %Load the first image
            LoadImg = double(imread(VarCell{1}.Analysis.FilePath{FilterNum,1}));
            
%             figure(94)
%             imagesc(LoadImg)
            
            %Initalize the Stack as uint16
            ImgSizeX = size(LoadImg,2);
            ImgSizeY = size(LoadImg,1);
            ImgStack = zeros(NumPos, size(LoadImg,1), size(LoadImg,2));
            %Filter image to remove pepper and salt noise
            LoadImg = medfilt2(LoadImg,'symmetric', [3 3 ]);
            
            %Other filtering processes
%             LoadImg = imopen(LoadImg,strel('disk',5));
%             
%             H = fspecial('gaussian',30,3);
%             LoadImg = imfilter(LoadImg,H,'symmetric');

            %Pass first image to stack
            ImgStack(1,:,:) = LoadImg;


            %Load all other images in the Stack
            for P = 2:NumPos
                LoadImg =double(imread(VarCell{P}{1}.Analysis.FilePath{FilterNum,1}));
                LoadImg = medfilt2(LoadImg,'symmetric', [3 3 ]);
                %Filter image to remove pepper and salt noise
                LoadImg = medfilt2(LoadImg,'symmetric', [3 3 ]);
                
                %Other filtering processes
                %             LoadImg = imopen(LoadImg,strel('disk',5));
                %
                %             H = fspecial('gaussian',30,3);
                %             LoadImg = imfilter(LoadImg,H,'symmetric');
                
                %Pass  image to stack
                ImgStack(P,:,:) = LoadImg;
            end
           
            %If there are more than 5 pos remove lower and 2 highest
            %intensities for each pixel
            if NumPos > 5
                %Sort the stack
                ImgStack = sort(ImgStack,1);
                %Calculate the mean Std and CV
                MeanStack = squeeze(mean(ImgStack(ceil(NumPos/10)+1:ceil(NumPos/2)+1,:,:),1));
                StdStack = squeeze(std(ImgStack(ceil(NumPos/10)+1:ceil(NumPos/2)+1,:,:),0,1));
                CVStack = StdStack./MeanStack;
            else
                %Calculate the mean Std and CV for all images in the stack
                MeanStack = squeeze(mean(ImgStack,1));
                StdStack = squeeze(std(ImgStack,0,1));
                CVStack = StdStack./MeanStack;
            end
%              figure(96)
%             imagesc(MeanStack)
%             
%             MidLine_Nothing = MeanStack(1000,:);
%             figure(98)
%             plot(MeanStack(1000,:))
            
            
%             %% Reinhard code
%             
%             if 0
%                 tic
%                 windowSize = 10;
%                 
%                 startx = [1:ImgSizeX] - windowSize;
%                 startx(startx<1) = 1;
%                 endx = [1:ImgSizeX] + windowSize;
%                 endx(endx>ImgSizeX) = ImgSizeX;
%                 
%                 starty = [1:ImgSizeY] - windowSize;
%                 starty(starty<1) = 1;
%                 endy = [1:ImgSizeY] + windowSize;
%                 endy(endy>ImgSizeY) = ImgSizeY;
%                 
%                 % smooth CVimage and detect threshold for allowed variation
%                 corrCVStack = zeros(ImgSizeY, ImgSizeX);
%                 corrFactor = zeros(ImgSizeY,ImgSizeX);
%                 for i=1:ImgSizeX
%                     for j=1:ImgSizeY
%                         pixels = CVStack(starty(j):endy(j),startx(i):endx(i));
%                         corrCVStack(j,i)=median(pixels(:));
%                     end
%                 end
%                 figure(201)
%                 imagesc(corrCVStack)
%                 
%                 cutoff = prctile (corrCVStack(:),70)
%                 
%                 % eliminate all pixels with higher variation than threshold, smooth
%                 % others
%                 
%                 for i=1:ImgSizeX
%                     for j=1:ImgSizeY
%                         before(j,i)= MeanStack(j,i);
%                         if (corrCVStack(j,i) < cutoff)
%                             pixels = MeanStack(starty(j):endy(j),startx(i):endx(i));
%                             corrFactor(j,i)=median(pixels(:));
%                         else
%                             corrFactor(j,i)=NaN;
%                         end
%                     end
%                 end
%                 figure(202)
%                 imagesc(corrFactor)
%                 
%                 % fill in gaps in smoothened image
%                 
%                 for  i=1:ImgSizeX
%                     for j=1:ImgSizeY
%                         if isnan(corrFactor(j,i))
%                             dx1=0;
%                             dx2=0;
%                             dy1=0;
%                             dy2=0;
%                             
%                             % try
%                             for k=1:ImgSizeX-i
%                                 if ~isnan(corrFactor(j,i+k))
%                                     dx1=k;
%                                     break
%                                 end
%                             end
%                             % end
%                             % try
%                             for k=1:i-1
%                                 if ~isnan(corrFactor(j, i-k))
%                                     dx2=k;
%                                     break
%                                 end
%                             end
%                             % end
%                             % try
%                             for k=1:ImgSizeY-j
%                                 if ~isnan(corrFactor(j+k,i))
%                                     dy1=k;
%                                     break
%                                 end
%                             end
%                             % end
%                             % try
%                             for k=1:j-1
%                                 if ~isnan(corrFactor(j-k,i))
%                                     dy2=k;
%                                     break
%                                 end
%                             end
%                             % end
%                             
%                             corrFactor(j,i)= nanmedian([corrFactor(j,i+dx1:endx(i+dx1)) corrFactor(j,startx(i-dx2):i-dx2) ...
%                                 corrFactor(j+dy1:endy(j+dy1),i)' corrFactor(starty(j-dy2):j-dy2,i)']);
%                             
%                             
%                         end
%                     end
%                 end
%                 
%                 figure(203)
%                 imagesc(corrFactor)
%                 
%                 %    final smoothening of the correction factor
%                 
%                 for i=1:ImgSizeX
%                     for j=1:ImgSizeY
%                         pixels = corrFactor(starty(j):endy(j),startx(i):endx(i));
%                         corrFactor(i,j)=median(pixels(:));
%                     end
%                 end
%                 
%                 figure(204)
%                 imagesc(corrFactor)
%                 
%                 toc
%                 drawnow
%             end
            %% my code

            

            %Median filtering of coefVariation image
            medCVStack = medfilt2(CVStack,'symmetric', [WindowSize WindowSize]);

%             figure(101)
%             imagesc(medCVStack)
%             
            
            %Median filter of MeanStack image
            medStack = medfilt2(MeanStack,'symmetric', [WindowSize WindowSize]);
            medStackFull = medStack;            
%                                     H = fspecial('gaussian',20,3);
%                     GaussStack = imfilter(medStack,H,'symmetric');
%                                   figure(122)
%             imagesc(GaussStack)
%             medStack = GaussStack;
            %Set 30% of pixels with a largest coef of vaiation to NaN
            cutoff = prctile(medCVStack(:),70);
            medStack(medCVStack>cutoff) = NaN;
%             figure(102)
%             imagesc(medStack)
            

            %Interpolate points with large CV
            %Create grid for the whole image
            [Xgrid,Ygrid] = meshgrid([1:ImgSizeY],[1:ImgSizeX]);
            
            %Find not NAN pixels
            Ind = find(~isnan(medStack));
            %Get value and XY coordinates for each pixel
            Zval = medStack(Ind);
            [Yind, Xind] = ind2sub(size(medStack),Ind);
            %Interpolate value of NAN pixels
            InterpStack = griddata(Xind, Yind, Zval,Xgrid,Ygrid,'linear');
            %Same as griddata maybe faster
            %             F = TriScatteredInterp(Xind, Yind, Zval);
            %             InterpStack = F(Xgrid,Ygrid);
%             figure(103)
%             imagesc(InterpStack)
            
            %Gaussian filter of Interpolated image
%             H = fspecial('gaussian',WindowSize,3);
%             GaussInterp = imfilter(InterpStack,H,'symmetric');
            % Median filter the interpolated image
            medInterp = medfilt2(InterpStack,'symmetric', [WindowSize WindowSize]);
%             figure(104)
%             imagesc(medInterp)
            %fill NAN values with med stack values
            Ind = find(isnan(medInterp));
            medInterp(Ind) = medStackFull(Ind);
            
            %Calculate flatness correction
            Flatness.(VarCell{1}.Analysis.LoadImgOut{I}) = medInterp./max(medInterp(:));
            
            %Save Flatness to file
            FlatFile = fullfile(fileparts(VarCell{1}.Analysis.OutPath), 'Flatness.mat');
            save(FlatFile, 'Flatness')
%             figure(105)
%             plot(Flatness(1000,:))
%             hold all
%             plot(Flatness(500,:))
%             plot(Flatness(1500,:))
            
            figure(100+I)
            subplot(2,2,1); imagesc(medStackFull); title('Median of image stack')
            subplot(2,2,2); imagesc(Flatness.(VarCell{1}.Analysis.LoadImgOut{I})); title('Flatness')
            subplot(2,2,4); plot(Flatness.(VarCell{1}.Analysis.LoadImgOut{I})(1000,:), 'r'); title('Linescan')
            hold all
            plot(Flatness.(VarCell{1}.Analysis.LoadImgOut{I})(500,:), 'g')
            plot(Flatness.(VarCell{1}.Analysis.LoadImgOut{I})(1500,:), 'b')
            axis tight
            drawnow
            
        end
        
        for P = 1:length(VarCell)
            VarCell{P}.Analysis.FlatFile = FlatFile;
        end
    end
end
