function [NuclImg, CellImg, GoodNucl] = ResolveNuclAssignment(NuclInitImg, NuclImg, CellImg,ProblemCell)
G = 0;
for P =1:length(ProblemCell)
    P = P
    NuclSize = [];
    Intersect_CN = [];
    IntersectSize = [];
    OverlapRatio = [];
    %Check if SmallOverlaps have been removed in NuclInit Image
    CellInd = find(CellImg == ProblemCell(P));
    
     correspondingNuclLabel = unique(NuclInitImg(CellInd));
        %Remove zeros pixels
        correspondingNuclLabel = correspondingNuclLabel(correspondingNuclLabel>0);
        if length(correspondingNuclLabel) == 1
             %Identify pixel indices
            NuclInd = find(NuclInitImg == correspondingNuclLabel);
            % Find overlap with cell indices
            OverlapInd = intersect(NuclInd,CellInd);
            %Set Overlap pixels to Nucleus Img to C
            NuclImg(OverlapInd) =ProblemCell(P);
            
            G = G+1;
            GoodNucl(G) = ProblemCell(P);
            continue
        end
        
        % Measure size of nucl and their overlap with the Cell
        for N = 1:length(correspondingNuclLabel)
            NuclInd = find(NuclInitImg == correspondingNuclLabel(N));
            NuclSize(N) = length(NuclInd);
            Intersect_CN{N} = (intersect(CellInd, NuclInd));
            IntersectSize(N) = length(Intersect_CN{N});
            OverlapRatio(N) = IntersectSize(N)/NuclSize(N);
        end
        
        %Test is only One nucleus is completely in the cell
        if length(find(OverlapRatio >= 0.9)) ==1
           GoodNucl(P) =  correspondingNuclLabel(OverlapRatio >= 0.9);
        elseif length(find(OverlapRatio >= 0.9))> 1
            [~, SizeSortInd] = sort(IntersectSize, 'descend');
            GoodNucl(P) = correspondingNuclLabel(SizeSortInd(1));
        else
            %This cell does not have a nucleus fully inside remove it
            CellImg(CellInd) = 0;
            continue
        end
        
           
         NuclInd = find(NuclInitImg == GoodNucl(P));
            % Find overlap with cell indices
            OverlapInd = intersect(NuclInd,CellInd);
            %Set Overlap pixels to Nucleus Img to C
            NuclImg(OverlapInd) = ProblemCell(P);
         G = G+1;
            GoodNucl(G) = ProblemCell(P);
end
