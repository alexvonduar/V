function [lsize] = line_size(seglines)
lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
end