function [pts]=createLogSection(start_Pts,end_Pts,NumSpaces,direction)
length=end_Pts-start_Pts;
x = linspace(0, 1, NumSpaces).^(1/3) * length;
pts=x+start_Pts;
if(strcmp(direction,'right'))
    x=-fliplr(x);
    pts=x+start_Pts+length;
end
end
