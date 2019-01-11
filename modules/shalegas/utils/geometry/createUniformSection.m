function [pts]=createUniformSection(start_Pts,end_Pts,NumSpaces)
length=end_Pts-start_Pts;
x = linspace(0, 1, NumSpaces)* length;
pts=x+start_Pts;
end