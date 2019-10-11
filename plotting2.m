analysis = load('analysis_result_c.txt');
% recurdyn = load('../recurdyn_result.txt');
recurdyn = load('../recurdyn_result_collision.txt');

figure
set(gcf,'Color',[1,1,1])
subplot(2,2,1)
plot(recurdyn(:,2), recurdyn(:,3),'LineWidth',1.5)
hold on
plot(analysis(:,1), analysis(:,2), '--', 'LineWidth',1.5)
grid on
xlabel('Time[sec]')
ylabel('Position[rad]')
title('q1 position')
set(gca,'FontSize',13)

subplot(2,2,2)
plot(recurdyn(:,2), recurdyn(:,4),'LineWidth',1.5)
hold on
plot(analysis(:,1), analysis(:,3), '--', 'LineWidth',1.5)
grid on
xlabel('Time[sec]')
ylabel('Position[meter]')
title('end x')
set(gca,'FontSize',13)
legend('RecurDyn','Analysis')

subplot(2,2,3)
plot(recurdyn(:,2), recurdyn(:,5),'LineWidth',1.5)
hold on
plot(analysis(:,1), analysis(:,4), '--', 'LineWidth',1.5)
grid on
xlabel('Time[sec]')
ylabel('Position[meter]')
title('end y')
set(gca,'FontSize',13)

subplot(2,2,4)
plot(recurdyn(:,2), recurdyn(:,6),'LineWidth',1.5)
hold on
plot(analysis(:,1), analysis(:,5), '--', 'LineWidth',1.5)
grid on
xlabel('Time[sec]')
ylabel('Position[meter]')
title('end z')
set(gca,'FontSize',13)

figure
set(gcf, 'Color', [1,1,1])
plot(analysis(:,1), analysis(:,6), 'LineWidth', 1.5)
grid on
xlabel('Time[sec]')
ylabel('Residual [Nm]')
set(gca,'FontSize',13)