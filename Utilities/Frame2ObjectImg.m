%close all
Var = Var_F1;

Object = 'Nucl';

NumObj = size(Var.Measurements.(Object).PixelList,2);
ObjImg = zeros(Var.Analysis.ImgSize);
Tpt = 4;


%Path = '/Volumes/DMF-1/GROUPS/gr_Pelet/Serge/2016/161007/1007_ySP643_ySP663_vs_ySP711_1/Pos2';
Path = '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2016/161007/1007_ySP643_ySP663_vs_ySP711_1/Pos0';
CFPname = num2str(Tpt-1, 'img_%09d_06-CFPtriple_000.tif')
CFPimg = double(imread(fullfile(Path, CFPname)));
CFPlim = [800, 3000];

CFPnorm = (CFPimg-CFPlim(1))./(CFPlim(2)-CFPlim(1));
CFPnorm(CFPnorm>1) = 1;
CFPnorm(CFPnorm<0) = 0;


for obj = 1:NumObj
    ObjImg(Var.Measurements.(Object).PixelList{Tpt, obj}) = obj;
end

figure(1)
imagesc(ObjImg)

ObjImgNorm = zeros(Var.Analysis.ImgSize);
ObjImgNorm(ObjImg>0) = 0.5;


RGBimg(:,:,3) = CFPnorm;
RGBimg(:,:,1) = ObjImgNorm;
figure(2)
image(RGBimg)

