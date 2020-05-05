ImgObj = zeros(Var.Analysis.ImgSize(1),Var.Analysis.ImgSize(2));

for O = 1:length(Var.Measurements.Nucl.PixelList)
    ImgObj(Var.Measurements.Nucl.PixelList{O}) = 1;
end

DAPI = double(imread('/Volumes/DMF/GROUPS/gr_Pelet/Min_B/2016/160119/MEK2_2NLS_mCherry_SC_2/Pos0/img_000000000_09-DAPIquad_000.tif'));
DAPI = (DAPI-150)./(6000-150);
DAPI(DAPI>1) = 1;
DAPI(DAPI<0) = 0;


RGB(:,:,1) = ImgObj;
RGB(:,:,2) = DAPI;
RGB(:,:,3) = 0;

figure(100)
image(RGB)