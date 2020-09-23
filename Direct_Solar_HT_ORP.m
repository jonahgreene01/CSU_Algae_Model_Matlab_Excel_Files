function [Q_Solar]= Direct_Solar_HT_ORP(GHI)
%This function defines the heat flux from direct solar radiation

f_a = 0.025; %fraction of sunlight converted to chemical energy during photosynthesis
Q_Solar = (1-f_a) * GHI*0.82; 

% Source: Yadala and Cremaschi, 2016

end

