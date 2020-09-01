function RecImage_selection = VirtualReal_ImageRec_selection (RecImage,  rect1, filterShape)
switch filterShape
    case 1
        [R,C]=size(Mag_RecImage1);
        line([round(C) round(C)],[0 round(R)])
        line([0 round(C)],[round(R) round(R)])
        h = imellipse(gca,[150,150,RR2,RR2]);
        vertices = wait(h);
        % load ('vertices.mat');
        XR=vertices(:,1);
        YR=vertices(:,2);
        Mask = poly2mask (XR,YR,R,C);
        Mag_RecImage_selection1=Mask.*Mag_RecImage1;
        Mag_RecImage_selection1 = Mag_RecImage_selection1(floor(min(YR)):floor(max(YR)),floor(min(XR)):floor(max(XR)));
        Mag_RecImage_selection2=Mask.*Mag_RecImage2;
        Mag_RecImage_selection2 = Mag_RecImage_selection2(floor(min(YR)):floor(max(YR)),floor(min(XR)):floor(max(XR)));
        
        Phase_RecImage1=Mask.*Phase_RecImage1;
        Phase_RecImage1 = Phase_RecImage1(floor(min(YR)):floor(max(YR)),floor(min(XR)):floor(max(XR)));
        Phase_RecImage2=Mask.*Phase_RecImage2;
        Phase_RecImage2 = Phase_RecImage2(floor(min(YR)):floor(max(YR)),floor(min(XR)):floor(max(XR)));
        
    case 2
        
        % load ('rect1.mat');
        RecImage_selection = RecImage(round(rect1(2)):round(rect1(2)+rect1(4)),round(rect1(1)):round(rect1(1)+rect1(3)));;
        Phase_RecImage = angle(RecImage_selection);
        Mag_RecImage_selection = abs(RecImage_selection);
        
end
end