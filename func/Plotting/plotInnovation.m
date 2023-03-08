function [] = plotInnovation(FilteredData)

figure('Position', [10, 10, 800, 800]);
sigma = 1;

% Build plot title
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Measurement innovations with \pm1\sigma bounds', ...
  'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', ...
  'Center', 'VerticalAlignment', 'Bottom' ) ;

clr_inn = [0.0 0.0 0.0]; % nearly black
clr_conf  = [0.7 0.4 0.4]; % grey bounds
lw_inn = 1.0; 
lw_conf = 2.0; % line width

%% plot innovation x-position
t = FilteredData{1}.Measurement.x.time;
innovation = FilteredData{1}.Measurement.x.innovation;
CovarianceBound = sigma*FilteredData{1}.Measurement.x.std;

subplot(4,4,1)
hold on
    grid on
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn, 'LineWidth', 1)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('x-position [m]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end
    
%% plot innovation y-position
t = FilteredData{1}.Measurement.y.time;
innovation = FilteredData{1}.Measurement.y.innovation;
CovarianceBound = sigma*FilteredData{1}.Measurement.y.std;

subplot(4,4,2)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('y-position [m]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end
    
%% plot innovation altitude
t = FilteredData{1}.Measurement.h.time;
innovation = FilteredData{1}.Measurement.h.innovation;
CovarianceBound = sigma*FilteredData{1}.Measurement.h.std;

subplot(4,4,3)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('altitude [m]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end   

%% plot innovation xdot
t = FilteredData{1}.Measurement.xdot.time;
innovation = FilteredData{1}.Measurement.xdot.innovation;
CovarianceBound = sigma*FilteredData{1}.Measurement.xdot.std;

subplot(4,4,4)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('xdot [m/s]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end
    
%% plot innovation ydot
t = FilteredData{1}.Measurement.ydot.time;
innovation = FilteredData{1}.Measurement.ydot.innovation;
CovarianceBound = sigma*FilteredData{1}.Measurement.ydot.std;

subplot(4,4,5)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('ydot [m/s]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end
    
    
%% plot innovation hdot
t = FilteredData{1}.Measurement.hdot.time;
innovation = FilteredData{1}.Measurement.hdot.innovation;
CovarianceBound = sigma*FilteredData{1}.Measurement.hdot.std;

subplot(4,4,6)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('hdot [m/s]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end    
    
%% plot innovation phi
t = FilteredData{1}.Measurement.phi.time;
innovation = FilteredData{1}.Measurement.phi.innovation*180/pi;
CovarianceBound = sigma*FilteredData{1}.Measurement.phi.std*180/pi;

subplot(4,4,7)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('phi [deg]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end
    
%% plot innovation theta
t = FilteredData{1}.Measurement.theta.time;
innovation = FilteredData{1}.Measurement.theta.innovation*180/pi;
CovarianceBound = sigma*FilteredData{1}.Measurement.theta.std*180/pi;

subplot(4,4,8)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('theta [deg]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end
    
%% plot innovation psi
t = FilteredData{1}.Measurement.psi.time;
innovation = FilteredData{1}.Measurement.psi.innovation*180/pi;
CovarianceBound = sigma*FilteredData{1}.Measurement.psi.std*180/pi;

subplot(4,4,9)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('psi [deg]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end
    
%% plot innovation Vtas
t = FilteredData{1}.Measurement.Vtas.time;
innovation = FilteredData{1}.Measurement.Vtas.innovation;
CovarianceBound = sigma*FilteredData{1}.Measurement.Vtas.std;

subplot(4,4,10)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('Vtas [m/s]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end   

%% plot innovation alpha
t = FilteredData{1}.Measurement.alpha.time;
innovation = FilteredData{1}.Measurement.alpha.innovation*180/pi;
CovarianceBound = sigma*FilteredData{1}.Measurement.alpha.std*180/pi;

subplot(4,4,11)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('alpha [deg]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end

%% plot innovation beta
t = FilteredData{1}.Measurement.beta.time;
innovation = FilteredData{1}.Measurement.beta.innovation*180/pi;
CovarianceBound = sigma*FilteredData{1}.Measurement.beta.std*180/pi;

subplot(4,4,12)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('beta [deg]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end  

%% plot innovation alpha_L
t = FilteredData{1}.Measurement.alpha_L.time;
innovation = FilteredData{1}.Measurement.alpha_L.innovation*180/pi;
CovarianceBound = sigma*FilteredData{1}.Measurement.alpha_L.std*180/pi;

subplot(4,4,13)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('alpha_L [deg]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end

%% plot innovation alpha_R
t = FilteredData{1}.Measurement.alpha_R.time;
innovation = FilteredData{1}.Measurement.alpha_R.innovation*180/pi;
CovarianceBound = sigma*FilteredData{1}.Measurement.alpha_R.std*180/pi;

subplot(4,4,14)
hold on     
    grid on 
    plot(t,  innovation, 'Color', clr_inn, 'LineWidth', lw_inn)
    plot(t,  CovarianceBound ,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    plot(t, -CovarianceBound,'--', 'Color', clr_conf, 'LineWidth', lw_conf)
    ylabel('alpha_R [deg]')
    xlabel('time [s]')
    if isempty(t) ~= 1
    ymin = min([min(innovation(500:end)) , min(-CovarianceBound(500:end))]);
    ymax = max([max(innovation(500:end)) , max(CovarianceBound(500:end))]);
    ylim([1.1*ymin 1.1*ymax])   
    end
    
    
%% Save figure
% saveas(gcf,'Innovation.png')
print('Innovation','-depsc')
end