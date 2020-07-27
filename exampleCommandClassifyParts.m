function exampleCommandClassifyParts(coordinator)
%
%CommandClassifyParts Classify the parts to determine where to place them
%   This command classifies the detected parts using a numeric type: type 1
%   or type 2.
%
% Copyright 2020 The MathWorks, Inc.

    % In this method, the classification is assumed to be known. In
    % practice, this may be replaced by more complex classification
    % tools. For example, the robot can use camera images or point
    % clouds  from 3D scanners to classify objects in the scene.
    
    %part1.color = 'w'  mixing beaker
    %part3.color = 'b' - distilled water
    %part2.color = 'r' -  NaCl
    %part4.color = 'g' - mixing platform
    %part5.color = 'y' - mixing platform
    
    % type 1 - placed on magnetic mixer
    % type 2 - poured over mixing beaker and then placed on the middle
    % shelf
    % type 3 - mixing beaker placed back in center of work table
    
    if ~isempty(coordinator.DetectedParts)
        coordinator.DetectedParts{3}.type = 1;
        coordinator.DetectedParts{2}.type = 2;
        coordinator.DetectedParts{1}.type = 3;
        coordinator.DetectedParts{1}.type = 4;
    end

   % Trigger Stateflow chart Event
   coordinator.FlowChart.partsClassified;       
end
