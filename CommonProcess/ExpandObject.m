function Var = ExpandObject(Var, CallNum)
tic
if nargin == 1
    CallNum = 1;
end

Debug =0;

%Assigne Ref and input images
ExpObjSmall = Var.Analysis.ExpObjSmall{CallNum};

ExpObjOUT = Var.Analysis.ExpObjOUT{CallNum};
Growth = Var.Analysis.ExpGrowth{CallNum};
Border = Var.Analysis.ExpBorder{CallNum};
Hole = Var.Analysis.ExpHole{CallNum};




%% Set Images for small and large Object

SmallObj =  Var.Img.(ExpObjSmall);

if isfield(Var.Analysis, 'ExpObjLarge') && (length(Var.Analysis.ExpObjLarge)  >= CallNum) && ~isempty(Var.Analysis.ExpObjLarge{CallNum})
    ExpObjLarge = Var.Analysis.ExpObjLarge{CallNum};
    LargeObj =  Var.Img.(ExpObjLarge);
else
    %Use Watershed to define large objects
    SmallBW = im2bw(SmallObj, 0.5);
    Dist = bwdist(SmallBW);
    LargeObj = watershed(Dist);
end

%Get Label from Small Obj
SmallLabel = Var.Measurements.(ExpObjSmall).SegLabel;
ObjProps = regionprops(SmallObj,'Centroid');

%% Grow object
ExpandObj = zeros(size(SmallObj));
for O = 1:max(SmallObj(:))
    %Create center obj and Around obj images containing only the objects of
    %interest
    CenterObj = zeros(size(SmallObj));
    AroundObj = zeros(size(SmallObj));
    %Place Small Object in CentreObj
    CenterObj(SmallObj == O) = 1;
    
    %Get LargeObj Label and place in AroundObj
    LargeObjNum = unique(LargeObj(SmallObj == O));
    AroundObj(LargeObj == LargeObjNum) = 1;
    %AroundObj(LargeObj == LargeObjNum(1)) = 1;
    
    %GrowCenterObj
    SE = strel('disk',Growth);
    GrownObj = imdilate(CenterObj, SE);

    %Make Hole in Objects
    if strcmpi(Hole, 'Yes')
        SE = strel('disk',Border);
        RemObj = imdilate(CenterObj, SE);
        GrownObj(RemObj == 1) = 0;
    end
    
    %Remove Pixel which are outside of Larger object
    GrownObj(AroundObj == 0) = 0;
    
    %Label objects
    [GrownObj, NumGrown] = bwlabel(GrownObj);
    if NumGrown == 0
        %% pb we don't have area surrounding the object
        %Set GrownObj to a single pixel at the center of the central
        %object
        PixList = find(CenterObj == 1);
        GrownObj(PixList(1)) = 1;
    end
%     SE = strel('disk',2);
%         GrownObj = imopen(GrownObj, SE);
%     
    
    %Add Object to Expand Object
    ExpandObj(GrownObj > 0) = O;
    Var.Measurements.(ExpObjOUT).CenterX(O,1) = ObjProps(O).Centroid(1);
    Var.Measurements.(ExpObjOUT).CenterY(O,1) = ObjProps(O).Centroid(2);
    Var.Measurements.(ExpObjOUT).PixelList{1,O} = find(GrownObj == 1);
    
end

Var.Measurements.(ExpObjOUT).SegLabel = SmallLabel;


%% Display %%%
if strcmp(Var.Figure.Display, 'on')
    FigNum = find(strcmp(Var.Figure.List, 'ExpObjSmall'));
    figure(FigNum(CallNum))
    subplot(2,2,1); imagesc(SmallObj);title('Small Object')
    subplot(2,2,2); imagesc(LargeObj); title('Large Object');
    subplot(2,2,3); imagesc(ExpandObj); title('Expanded Object Object');
    
    Obj = zeros(size(SmallObj));
    Obj(ExpandObj>0) = 1;
    BothObjRGB(:,:,1) = Obj;
    Obj = zeros(size(SmallObj));
    %Obj(LargeObj>0) = 1;
    BothObjRGB(:,:,2) = Obj;
    Obj = zeros(size(SmallObj));
    Obj(SmallObj>0) = 1;
    BothObjRGB(:,:,3) = Obj;
    subplot(2,2,4); image(BothObjRGB); title('Overlay');
    
end


%% save %%%
Var.Img.(ExpObjOUT) = ExpandObj;


%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;

