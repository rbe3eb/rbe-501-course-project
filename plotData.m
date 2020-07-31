function plotData(T,X,Xdes,eePos,eePosDes,method,space)
    % Joint Position
    figure(2)
    for n=1:7
        plot(T, X(:,n));
        hold on
    end
    xlim([0 T(end)])
    legend('q_1', 'q_2', 'q_3')
    xlabel('time [sec]')
    ylabel('position [mm]')
    title(['Actual Joint Position', ' (', method, ')', newline, space])
    hold off
    
    % Joint Velocity
    figure(3)
    for n=8:14
        plot(T, X(:,n));
        hold on
    end
    xlim([0 T(end)])
    legend('dq_1', 'dq_2', 'dq_3')
    xlabel('time [sec]')
    ylabel('velocity [mm/sec]')
    title(['Actual Joint Velocity', ' (', method, ')', newline, space])
    hold off

    % Joint Position
    figure(4)
    for n=1:7
        plot(T, Xdes(:,n));
        hold on
    end
    xlim([0 T(end)])
    legend('q_1', 'q_2', 'q_3')
    xlabel('time [sec]')
    ylabel('position [mm]')
    title(['Desired Joint Position', ' (', method, ')', newline, space])
    hold off
    
    % Joint Velocity
    figure(5)
    for n=8:14
        plot(T, Xdes(:,n));
        hold on
    end
    xlim([0 T(end)])
    legend('dq_1', 'dq_2', 'dq_3')
    xlabel('time [sec]')
    ylabel('velocity [mm/sec]')
    title(['Desired Joint Velocity', ' (', method, ')', newline, space])
    hold off
    drawnow
    
    
    % End-Effector Position
    % x
    figure(6)
    subplot(1,3,1)
    plot(T,eePos(:,1))
    hold on
    plot(T,eePosDes(:,1))
    hold off
    xlabel('time [sec]')
    ylabel('x [mm]')
    legend(['Actual Motion', 'Desired Motion'])
    xlim([0 T(end)]);
    title('(x-axis)')
    %y
    subplot(1,3,2)
    plot(T,eePos(:,2))
    hold on
    plot(T,eePosDes(:,1))
    hold off
    xlabel('time [sec]')
    ylabel('y [mm]')
    legend(['Actual Motion', 'Desired Motion'])
    xlim([0 T(end)]);
    title('(y-axis)')
    %z
    subplot(1,3,3)
    plot(T,eePos(:,3))
    hold on
    plot(T,eePosDes(:,1))
    xlabel('time [sec]')
    ylabel('z [mm]')
    legend(['Actual Motion', 'Desired Motion'])
    xlim([0 T(end)]);
    title('(z-axis)')
    hold off
    sgtitle(['End-Effector Motion', ' (', method, ')', newline, space])
    
    
end