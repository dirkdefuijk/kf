function [] = plotInput(FilteredData)
%
figure('Position', [50, 50, 800, 800]);

% Build plot title
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Inputs (grey: raw, dark blue: reconstructed)', ...
  'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', ...
  'Center', 'VerticalAlignment', 'Bottom' ) ;

clr_raw            = [0.6 0.6 0.6]; % light grey
clr_reconstructed  = [0.1 0.1 0.6]; % blue plot line

%% plot acceleration in x-direction
t = FilteredData{1}.Input.u.ax.time;
raw = FilteredData{1}.Input.u.ax.raw;
filtered = FilteredData{1}.Input.u.ax.filtered;

subplot(3,3,1)
hold on
    grid on
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  filtered,'Color', clr_reconstructed)
    ylabel('A_x [m/s^2]')
    xlabel('time [s]')
    
%% plot acceleration in y-direction
t = FilteredData{1}.Input.u.ay.time;
raw = FilteredData{1}.Input.u.ay.raw;
filtered = FilteredData{1}.Input.u.ay.filtered;

subplot(3,3,2)
hold on
    grid on
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  filtered,'Color', clr_reconstructed)
    ylabel('A_y [m/s^2]')
    xlabel('time [s]')    
    
%% plot acceleration in z-direction
t = FilteredData{1}.Input.u.az.time;
raw = FilteredData{1}.Input.u.az.raw;
filtered = FilteredData{1}.Input.u.az.filtered;

subplot(3,3,3)
hold on
    grid on
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  filtered,'Color', clr_reconstructed)
    ylabel('A_z [m/s^2]')
    xlabel('time [s]')      
%% plot rotational rate p
t = FilteredData{1}.Input.u.p.time;
raw = FilteredData{1}.Input.u.p.raw;
filtered = FilteredData{1}.Input.u.p.filtered;

subplot(3,3,4)
hold on
    grid on
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  filtered,'Color', clr_reconstructed)
    ylabel('p [rad/s]')
    xlabel('time [s]') 
    
%% plot rotational rate q
t = FilteredData{1}.Input.u.q.time;
raw = FilteredData{1}.Input.u.q.raw;
filtered = FilteredData{1}.Input.u.q.filtered;

subplot(3,3,5)
hold on
    grid on
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  filtered,'Color', clr_reconstructed)
    ylabel('q [rad/s]')
    xlabel('time [s]')      
    
%% plot rotational rate r
t = FilteredData{1}.Input.u.r.time;
raw = FilteredData{1}.Input.u.r.raw;
filtered = FilteredData{1}.Input.u.r.filtered;

subplot(3,3,6)
hold on
    grid on
    plot(t,  raw, 'Color', clr_raw)
    plot(t,  filtered,'Color', clr_reconstructed)
    ylabel('r [rad/s]')
    xlabel('time [s]')
    
%% plot elevator deflection
t = FilteredData{1}.Input.delta.de.time;
raw = FilteredData{1}.Input.delta.de.raw;
filtered = FilteredData{1}.Input.delta.de.filtered;

subplot(3,3,7)
hold on
    grid on
    plot(t,  raw*180/pi, 'Color', clr_raw)
    if ~isempty(filtered)
        plot(t,  filtered*180/pi,'Color', clr_reconstructed)
    end   
    ylabel('\delta_e [deg]')
    xlabel('time [s]')
    
%% plot aileron deflection
t = FilteredData{1}.Input.delta.da.time;
raw = FilteredData{1}.Input.delta.da.raw;
filtered = FilteredData{1}.Input.delta.da.filtered;

subplot(3,3,8)
hold on
    grid on
    plot(t,  raw*180/pi, 'Color', clr_raw)
    if ~isempty(filtered)
        plot(t,  filtered*180/pi,'Color', clr_reconstructed)
    end   
    ylabel('\delta_a [deg]')
    xlabel('time [s]')
    
%% plot rudder deflection
t = FilteredData{1}.Input.delta.dr.time;
raw = FilteredData{1}.Input.delta.dr.raw;
filtered = FilteredData{1}.Input.delta.dr.filtered;

subplot(3,3,9)
hold on
    grid on
    plot(t,  raw*180/pi, 'Color', clr_raw)
    if ~isempty(filtered)
        plot(t,  filtered*180/pi,'Color', clr_reconstructed)
    end   
    ylabel('\delta_r [deg]')
    xlabel('time [s]')       
    
%% Save figure
print('Input','-depsc')
end