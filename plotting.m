function plotting

    global flag sim_flag

    if flag == 1
        recurdyn = load('ear3.txt');
        recurdyn_q = [recurdyn(:,2),recurdyn(:,3:8)];
        recurdyn_end = [recurdyn(:,2),recurdyn(:,9:14)];
        c_q = load('hj_kinematics_result.txt');
    elseif flag == 2
        recurdyn = load('ear3_motion_pick.txt');
        recurdyn_q = [recurdyn(:,2),recurdyn(:,3:8)];
        recurdyn_end = [recurdyn(:,2),recurdyn(:,9:14)];
        c_q = load('hj_inverse_kinematics_result.txt');
%     else
%         if sim_flag ~= 6
%             recurdyn_q = load(sprintf('inverse_kinematics_output_q%d.txt', sim_flag));
%             recurdyn_end = load(sprintf('inverse_kinematics_input_end%d.txt', sim_flag));
%         else
%             recurdyn_q = load('inverse_kinematics_output_q_pick.txt');
%             recurdyn_end = load('inverse_kinematics_input_end_pick.txt');
%         end
%         c_q = load('hj_kinematic_result.txt');
    end
%     recurdyn_q(1,:) = [];
%     recurdyn_q(end,:) = [];
%     recurdyn_end(1,:) = [];
%     recurdyn_end(end,:) = [];

    

    q1_sp = [c_q(:,1), c_q(:,2)];
    q2_sp = [c_q(:,1), c_q(:,3)];
    q3_sp = [c_q(:,1), c_q(:,4)];
    q4_sp = [c_q(:,1), c_q(:,5)];
    q5_sp = [c_q(:,1), c_q(:,6)];
    q6_sp = [c_q(:,1), c_q(:,7)];
    
    dlmwrite('q1_sp.txt',q1_sp,'\t');
    dlmwrite('q2_sp.txt',q2_sp,'\t');
    dlmwrite('q3_sp.txt',q3_sp,'\t');
    dlmwrite('q4_sp.txt',q4_sp,'\t');
    dlmwrite('q5_sp.txt',q5_sp,'\t');
    dlmwrite('q6_sp.txt',q6_sp,'\t');

    title_end = {'end point x','end point y','end point z','end point roll','end point pitch','end point yaw'};

    %% result compare plot
    figure
    set(gcf,'Color',[1,1,1])
    for i = 1 : 12
        subplot(4,3,i)
        if i <= 6
            plot(recurdyn_q(:,1), recurdyn_q(:,i+1),'b','LineWidth',2);
            hold on
            plot(c_q(:,1), c_q(:,i+1),'r--','LineWidth',2);
            grid on
            title(sprintf('q %d',i))
            xlabel('Time [s]')
            ylabel('Angle [deg]')
            set(gca,'FontSize',15)
        else
            if (i >= 7 && i <= 9)
                plot(recurdyn_end(:,1), recurdyn_end(:,i-5),'b','LineWidth',2);
                hold on
                plot(c_q(:,1), c_q(:,i+1),'r--','LineWidth',2);
            else
                plot(recurdyn_end(:,1), recurdyn_end(:,i-5)*180/pi,'b','LineWidth',2);
                hold on
                plot(c_q(:,1), c_q(:,i+1)*180/pi,'r--','LineWidth',2);
            end
            grid on
            title(title_end{i-6})
%             if  i == 10
%                 ylim([89 91])
%             end
            xlabel('Time [s]')
            if (i >= 7 && i <= 9)
                ylabel('Position [meter]')
            else
                ylabel('Angle [deg]')
            end
            set(gca,'FontSize',15)
        end
        if i == 3
            legend('RecurDyn','Analysis')
        end
    end

    %% error pot
%     figure
%     set(gcf,'Color',[1,1,1])
%     for i = 1 : 12
%         subplot(4,3,i)
%         if i <= 6
%             plot(recurdyn_q(:,1), recurdyn_q(:,i+1)*180/pi - c_q(:,i+1)*180/pi,'b','LineWidth',2.5);
%             grid on
%             title(sprintf('q %d',i))
%             xlabel('Time [s]')
%             ylabel('Angle [deg]')
%             set(gca,'FontSize',15)
%         else
%             if i >= 7 && i <= 9
%                 plot(recurdyn_end(:,1), recurdyn_end(:,i-5) - c_q(:,i+1),'b','LineWidth',2.5);
%             else
%                 plot(recurdyn_end(:,1), recurdyn_end(:,i-5)*180/pi - c_q(:,i+1)*180/pi,'b','LineWidth',2.5);
%             end
%             grid on
%             title(title_end{i-6})
%             xlabel('Time [s]')
%             if (i >= 7 && i <= 9)
%                 ylabel('Position [meter]')
%             else
%                 ylabel('Angle [deg]')
%             end
%             set(gca,'FontSize',15)
%         end
%     end

end