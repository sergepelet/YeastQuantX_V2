function Var = Img2MovieRGB(Var, CallNum)
tic
if nargin == 1
    CallNum = 1;
end
%Get image from var
Img = double(Var.Img.(Var.Analysis.MovieImg{CallNum}));
Var.Analysis.SaveAvi = 'No';
CurrentFrame = Var.Analysis.CurrentFrame;
if length(size(Img)) == 2
    %Normalize image between 0 and 1
    MaxInt = max(Img(:))*0.95;
    MinInt = min(Img(:))*1.05;
    NormImg = (Img-MinInt)./(MaxInt-MinInt);
    NormImg(NormImg <0) = 0;
    NormImg(NormImg >1) = 1;
end

%Create RGB image from Img based on defined colormap
RGBImg = zeros(size(Img,1),size(Img,2), 3);
if strcmpi(Var.Analysis.ColorMap{CallNum}, 'g')
    RGBImg(:,:,2) = NormImg;
elseif strcmpi(Var.Analysis.ColorMap{CallNum}, 'r')
    RGBImg(:,:,1) = NormImg;
elseif strcmpi(Var.Analysis.ColorMap{CallNum}, 'b')
    RGBImg(:,:,3) = NormImg;
elseif strcmpi(Var.Analysis.ColorMap{CallNum}, 'w')
    RGBImg(:,:,1) = NormImg;
    RGBImg(:,:,2) = NormImg;
    RGBImg(:,:,3) = NormImg;
elseif strcmpi(Var.Analysis.ColorMap{CallNum}, 'rgb')
    RGBImg = Img;
end
% figure(100)
% image(RGBImg)

%Check if objects have to be displayed on image
if isfield(Var.Analysis, 'MovieObjectDisplay') ...
        && length(Var.Analysis.MovieObjectDisplay) >= CallNum ...
        && ~isempty(Var.Analysis.MovieObjectDisplay{CallNum})...
        && isfield(Var.Measurements.(Var.Analysis.MovieObjectDisplay{CallNum}), 'PixelList')
    
    %Get object name and properties
    ObjectName = Var.Analysis.MovieObjectDisplay{CallNum};
    Measurements = Var.Measurements.(ObjectName);
    NbObjects = size(Measurements.PixelList,2);
    Label = Var.Measurements.(ObjectName).SegLabel;
    
    ObjImg = Var.Img.(ObjectName);
    ObjectsProps = regionprops(ObjImg,'Centroid');
    NumObjects = size(ObjectsProps,1);
    ObjBW = zeros(size(ObjImg));
    ObjBW(ObjImg >0) = 1;
    SatImg = ones(size(ObjImg));
    %ValueImg = bwperim(ObjBW);
    SE = strel('disk',2);
    ValueImg = ObjBW-imerode(ObjBW,SE);
    HueImg = zeros(size(ObjImg));
    NumColors = 10;
    for O = 1:NumObjects
        HueImg(ObjImg == O) = rem(Measurements.SegLabel(O), NumColors)/NumColors;
        ObjectLab_Pos(O,1) = Label(O);
        ObjectLab_Pos(O,2:3) = ObjectsProps(O).Centroid;
    end
    HSVImg(:,:,1) = HueImg;
    HSVImg(:,:,2) =SatImg;
    HSVImg(:,:,3) =ValueImg;
    RGBLabel = hsv2rgb(HSVImg);
    %     figure(101)
    % image(RGBLabel)
    %combine intensity and object images
    NotPerim =  imcomplement(ValueImg);
    NotPerim = NotPerim(:,:,[1,1,1]);
    RGBImg = RGBImg.*NotPerim + RGBLabel;
end

% figure(102)
% image(RGBImg)

%Set filename for movie
BaseFileName = Var.Analysis.OutPath;
if strcmp(Var.Analysis.SaveAvi, 'Yes')
    mov_name = [BaseFileName,'_',Var.Analysis.MovieImg{CallNum}, '.avi'];
    %Create avi object if first frame
    if strcmp(Var.Analysis.FrameAnalyzed,'First')
        %Make sure the avi object is not open for writing.
        if isfield(Var.Analysis, 'MovieObj')
            if length(Var.Analysis.MovieObj) >= CallNum
                if strcmp(get(Var.Analysis.MovieObj{CallNum}, 'CurrentState'), 'Open')
                    Var.Analysis.MovieObj{CallNum} = close(Var.Analysis.MovieObj{CallNum});
                end
            end
        end
        % Create avi object
        Var.Analysis.MovieObj{CallNum} = avifile(mov_name, 'fps', 2, 'Quality', 100, 'Compression', 'None');
    end
else
    % if strcmp(Var.Analysis.FrameAnalyzed,'First')
    mkdir([BaseFileName, '_', Var.Analysis.MovieImg{CallNum}]);
    %  end
    % Save each individual image in folder
    img_name = [BaseFileName,'_',Var.Analysis.MovieImg{CallNum},filesep,'fr', num2str(Var.Analysis.CurrentFrame,'%03d.jpeg')];
end

%Generate timestamp
if isfield(Var.Analysis, 'TimeStamp')
    FrameTime = Var.Analysis.TimeStamp- Var.Analysis.TimeZero;
else
    if strcmpi(Var.Experiment.TimeUnit, 'sec')
        TimeInterval = Var.Experiment.TimeStep*1000;
    elseif strcmpi(Var.Experiment.TimeUnit, 'min')
        TimeInterval = Var.Experiment.TimeStep*60*1000;
    end
    FrameTime = Var.Analysis.CurrentFrame*TimeInterval - Var.Experiment.TimeZero*TimeInterval;
end
%Transform time in to string
TimeStamp = gen_Tstamp(FrameTime);
TxtLoc = 1/10*size(Img);

%% Display Image from combined objects
if strcmp(Var.Figure.Display, 'on')
    FigNum = find((strcmp(Var.Figure.List, 'MakeMovie')));
    figure(FigNum(CallNum))
    % on PC set large image
    %if strcmp(filesep, '\')
    pos = [180; 251; 896; 672];
    set(FigNum(CallNum), 'Position', pos)
    %end
    %Set backgroung color to black
    set(FigNum(CallNum),'Color',[0 0 0])
    %set(FigNum, 'renderer', 'OpenGL')
    
    % Display image
    
    image(RGBImg)
    %Set axis positions such that the image fills the whole window
    set(gca,'Position',[0 0 1 1 ])
    
    %Display text on figure
    text(TxtLoc(1), TxtLoc(2),TimeStamp, 'Color', 'w')
    if isfield(Var.Analysis, 'MovieLabel') && strcmp(Var.Analysis.MovieLabel{CallNum}, 'Yes')
        for O = 1:size(ObjectLab_Pos,1)
            text(ObjectLab_Pos(O,2),ObjectLab_Pos(O,3), num2str(ObjectLab_Pos(O,1)), 'Color', 'w')
        end
    end
    set(gca,'ytick',[], 'xtick', []);
    drawnow
    if strcmp(Var.Analysis.SaveAvi, 'Yes')
        Var.Analysis.MovieObj{CallNum} = addframe(Var.Analysis.MovieObj{CallNum},FigNum(CallNum));
    else
        print('-djpeg',FigNum(CallNum),img_name);
    end
    
    
end
%pause(0.2)

% if strcmp(Var.Analysis.FrameAnalyzed,'Last') && strcmp(Var.Analysis.SaveAvi, 'Yes')
%     Var.Analysis.MovieObj{CallNum} = close(Var.Analysis.MovieObj{CallNum});
% end
%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function timestamp = gen_Tstamp(timediff);

Hrs = floor(abs(timediff)./(60*60*1000));
min_tot = rem(abs(timediff),60*60*1000);
Min = floor(min_tot./(60*1000));
sec_tot = rem(min_tot, 60*1000);
Sec = floor(sec_tot./(1000));
Msec = rem(sec_tot, 1000);

for i = 1:length(timediff)
    if timediff(i) >= 0
        sgn = ' ';
    else
        sgn = '-';
    end
    if max(Hrs) == 0
        timestamp(i,:) = [sgn, num2str(Min(i),'%02d'),''' :', num2str(Sec(i),'%02d'), ''''' :', num2str(Msec(i),'%03d')];
    else
        timestamp(i,:) = [sgn, num2str(Hrs(i),'%02d'),':',num2str(Min(i),'%02d'),''':', num2str(Sec(i),'%02d'), ''''''];
    end
end