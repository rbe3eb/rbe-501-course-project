function exampleCommandComputeGraspPose(coordinator) 

        coordinator.GraspPose = trvec2tform(coordinator.Parts{coordinator.NextPart}.centerPoint + [0 0 0.05])*axang2tform([0 1 0 pi]);

        % Trigger Stateflow chart Event
        %coordinator.FlowChart.nextAction; 
end