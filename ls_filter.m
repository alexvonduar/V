function seglines = ls_filter(thres_aligned, length_t, lsd)
seglines = lsd;
L0 = line_size(seglines);
seglinesC(1:size(seglines,1),1) = (seglines(:,3)+seglines(:,1))/2;
seglinesC(1:size(seglines,1),2) = (seglines(:,4)+seglines(:,2))/2;
seglinesJ(1:size(seglines,1),1) = -(seglines(:,4)-seglines(:,2))./L0;
seglinesJ(1:size(seglines,1),2) = (seglines(:,3)-seglines(:,1))./L0;
[L0, I] = sort(L0, 'descend');
seglines =  seglines(I,:);
seglinesC = seglinesC(I,:);
seglinesJ = seglinesJ(I,:);
i = 1;
while i < size(seglines,1)
    C = seglinesC(i,1:2);
    J = seglinesJ(i,:);
    X1 = abs((seglines((i+1):end,1:2)-repmat(C,size(seglines,1)-i,1))*J');
    X2 = abs((seglines((i+1):end,3:4)-repmat(C,size(seglines,1)-i,1))*J');
    I = intersect(find(X1 < thres_aligned),find(X2 < thres_aligned))+i;
    seglines(I,:) = [];
    L0(I) = [];
    seglinesC(I,:) = [];
    seglinesJ(I,:) = [];
    i = i+1;
end
i = 1;
while i <= size(seglines,1)
    if L0(i) < length_t/2
        seglines(i,:) = [];
        L0(i) = [];
        seglinesJ(i,:) = [];
        i = i-1;
    end
    i=i+1;
end
end