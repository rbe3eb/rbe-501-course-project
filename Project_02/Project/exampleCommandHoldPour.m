function exampleCommandHoldPour(coordinator)% This class is for internal use and may be removed in a future release


    if coordinator.NextPart.type == 1||2
               coordinator.NextPart.HoldConfig = true;
    else
        coordinator.DetectedParts{coordinator.NextPart}.HoldConfig = false;
    end
    