function plotData(T,X,method,space)
    % Joint Position
    figure()
    for n=1:3
        plot(T, X(:,n));
        hold on
    end
    xlim([0 T(end)])
    legend('q_1', 'q_2', 'q_3')
    xlabel('time [sec]')
    ylabel('position [mm]')
    title(['Joint Position', ' (', method, ')', newline, space])
    hold off
    
    % Joint Velocity
    figure()
    for n=4:6
        plot(T, X(:,n));
        hold on
    end
    xlim([0 T(end)])
    legend('dq_1', 'dq_2', 'dq_3')
    xlabel('time [sec]')
    ylabel('velocity [mm/sec]')
    title(['Joint Velocity', ' (', method, ')', newline, space])
    hold off
    
%     % End-Effector Position
%     % x
%     figure()
%     subplot(1,3,1)
%     plot(T,eePos(:,1))
%     xlabel('time [sec]')
%     ylabel('x [mm]')
%     xlim([0 tme(end)]);
%     title('(x-axis)')
%     %y
%     subplot(1,3,2)
%     plot(T,eePos(:,2))
%     xlabel('time [sec]')
%     ylabel('y [mm]')
%     xlim([0 tme(end)]);
%     title('(y-axis)')
%     %z
%     subplot(1,3,3)
%     plot(T,eePos(:,3))
%     xlabel('time [sec]')
%     ylabel('z [mm]')
%     xlim([0 tme(end)]);
%     title('(z-axis)')
%     sgtitle(['End-Effector Position', ' (', method, ')', newline, space])
%     
%     
    drawnow
end