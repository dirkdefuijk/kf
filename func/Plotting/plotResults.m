function [] = plotResults(FilteredData)
%
figure('Position', [20, 20, 800, 800]);

% Build plot title
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Measurements (grey: raw, dark blue: reconstructed)', ...
  'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', ...
  'Center', 'VerticalAlignment', 'Bottom' ) ;

clr_raw            = [0.6 0.6 0.6]; % light grey
clr_reconstructed  = [0.1 0.1 0.6]; % blue plot line

%% plot x-position
t = FilteredData{1}.Measurement.x.time;
raw = FilteredData{1}.Measurement.x.raw;
reconstructed = FilteredData{1}.Measurement.x.reconstructed;

subplot(4,4,1)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('x-position [m]')
    xlabel('time [s]')

%% plot y-position
t = FilteredData{1}.Measurement.y.time;
raw = FilteredData{1}.Measurement.y.raw;
reconstructed = FilteredData{1}.Measurement.y.reconstructed;

subplot(4,4,2)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('y-position [m]')
    xlabel('time [s]')
    
%% plot altitude
t = FilteredData{1}.Measurement.h.time;
raw = FilteredData{1}.Measurement.h.raw;
reconstructed = FilteredData{1}.Measurement.h.reconstructed;

subplot(4,4,3)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('altitude [m]')
    xlabel('time [s]')
   
%% plot xdot
t = FilteredData{1}.Measurement.xdot.time;
raw = FilteredData{1}.Measurement.xdot.raw;
reconstructed = FilteredData{1}.Measurement.xdot.reconstructed;

subplot(4,4,4)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('xdot [m/s]')
    xlabel('time [s]')
    
%% plot ydot
t = FilteredData{1}.Measurement.ydot.time;
raw = FilteredData{1}.Measurement.ydot.raw;
reconstructed = FilteredData{1}.Measurement.ydot.reconstructed;

subplot(4,4,5)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('ydot [m/s]')
    xlabel('time [s]')
    
%% plot hdot
t = FilteredData{1}.Measurement.hdot.time;
raw = FilteredData{1}.Measurement.hdot.raw;
reconstructed = FilteredData{1}.Measurement.hdot.reconstructed;

subplot(4,4,6)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('hdot [m/s]')
    xlabel('time [s]')  
    
%% plot phi
t = FilteredData{1}.Measurement.phi.time;
raw = FilteredData{1}.Measurement.phi.raw*180/pi;
reconstructed = FilteredData{1}.Measurement.phi.reconstructed*180/pi;

subplot(4,4,7)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('phi [deg]')
    xlabel('time [s]')
    
%% plot theta
t = FilteredData{1}.Measurement.theta.time;
raw = FilteredData{1}.Measurement.theta.raw*180/pi;
reconstructed = FilteredData{1}.Measurement.theta.reconstructed*180/pi;

subplot(4,4,8)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('theta [deg]')
    xlabel('time [s]')    
    
%% plot psi
t = FilteredData{1}.Measurement.psi.time;
raw = FilteredData{1}.Measurement.psi.raw*180/pi;
reconstructed = FilteredData{1}.Measurement.psi.reconstructed*180/pi;

subplot(4,4,9)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('psi [deg]')
    xlabel('time [s]') 
    
%% plot Vtas
t = FilteredData{1}.Measurement.Vtas.time;
raw = FilteredData{1}.Measurement.Vtas.raw;
reconstructed = FilteredData{1}.Measurement.Vtas.reconstructed;

subplot(4,4,10)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('Vtas [m/s]')
    xlabel('time [s]')
    
%% plot alpha
t = FilteredData{1}.Measurement.alpha.time;
raw = FilteredData{1}.Measurement.alpha.raw*180/pi;
reconstructed = FilteredData{1}.Measurement.alpha.reconstructed*180/pi;

subplot(4,4,11)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('alpha [deg]')
    xlabel('time [s]')   
    
%% plot beta
t = FilteredData{1}.Measurement.beta.time;
raw = FilteredData{1}.Measurement.beta.raw*180/pi;
reconstructed = FilteredData{1}.Measurement.beta.reconstructed*180/pi;

subplot(4,4,12)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('beta [deg]')
    xlabel('time [s]') 

%% plot alpha_L
t = FilteredData{1}.Measurement.alpha_L.time;
raw = FilteredData{1}.Measurement.alpha_L.raw*180/pi;
reconstructed = FilteredData{1}.Measurement.alpha_L.reconstructed*180/pi;

subplot(4,4,13)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('alpha_L [deg]')
    xlabel('time [s]')   

%% plot alpha_R
t = FilteredData{1}.Measurement.alpha_R.time;
raw = FilteredData{1}.Measurement.alpha_R.raw*180/pi;
reconstructed = FilteredData{1}.Measurement.alpha_R.reconstructed*180/pi;

subplot(4,4,14)
hold on
    grid on 
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  reconstructed ,'Color', clr_reconstructed)
    ylabel('alpha_R [deg]')
    xlabel('time [s]')   


%% Save figure
print('Measurements','-depsc')
end