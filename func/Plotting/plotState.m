function [] = plotState(FilteredData)
%
figure('Position', [30, 30, 800, 800]);
sigma = 1;

% Build plot title
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'State 1/3 (dark blue: state, grey dash: \pm1\sigma confidence bound)', ...
  'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', ...
  'Center', 'VerticalAlignment', 'Bottom' ) ;

clr_est = [0.1 0.1 0.6];
clr_conf = [0.5 0.5 0.5];

%% plot x-position
t = FilteredData{1}.State.x.time;
value = FilteredData{1}.State.x.data;
bounds = FilteredData{1}.State.x.std;

subplot(4,3,1)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('x-position [m]')
    xlabel('time [s]')
    
%% plot y-position
t = FilteredData{1}.State.y.time;
value = FilteredData{1}.State.y.data;
bounds = FilteredData{1}.State.y.std;

subplot(4,3,2)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('y-position [m]')
    xlabel('time [s]')  
    
%% plot z-position
t = FilteredData{1}.State.z.time;
value = FilteredData{1}.State.z.data;
bounds = FilteredData{1}.State.z.std;

subplot(4,3,3)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('z-position [m]')
    xlabel('time [s]')    
         
%% plot body velocity u
t = FilteredData{1}.State.u.time;
value = FilteredData{1}.State.u.data;
bounds = FilteredData{1}.State.u.std;

subplot(4,3,4)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('u [m/s]')
    xlabel('time [s]')
    
%% plot body velocity v
t = FilteredData{1}.State.v.time;
value = FilteredData{1}.State.v.data;
bounds = FilteredData{1}.State.v.std;

subplot(4,3,5)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('v [m/s]')
    xlabel('time [s]')
    
%% plot body velocity w
t = FilteredData{1}.State.w.time;
value = FilteredData{1}.State.w.data;
bounds = FilteredData{1}.State.w.std;

subplot(4,3,6)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('w [m/s]')
    xlabel('time [s]')
    
%% plot euler angle phi
t = FilteredData{1}.State.phi.time;
value = FilteredData{1}.State.phi.data*180/pi;
bounds = FilteredData{1}.State.phi.std*180/pi;

subplot(4,3,7)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('phi [deg]')
    xlabel('time [s]')
    
%% plot euler angle theta
t = FilteredData{1}.State.theta.time;
value = FilteredData{1}.State.theta.data*180/pi;
bounds = FilteredData{1}.State.theta.std*180/pi;

subplot(4,3,8)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('theta [deg]')
    xlabel('time [s]')  
    
%% plot euler angle psi
t = FilteredData{1}.State.psi.time;
value = FilteredData{1}.State.psi.data*180/pi;
bounds = FilteredData{1}.State.psi.std*180/pi;

subplot(4,3,9)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('psi [deg]')
    xlabel('time [s]')       

%% plot wind velocity Wxe
t = FilteredData{1}.State.Wxe.time;
value = FilteredData{1}.State.Wxe.data;
bounds = FilteredData{1}.State.Wxe.std;

subplot(4,3,10)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('Wxe [m/s]')
    xlabel('time [s]')
     
%% plot wind velocity Wye
t = FilteredData{1}.State.Wye.time;
value = FilteredData{1}.State.Wye.data;
bounds = FilteredData{1}.State.Wye.std;

subplot(4,3,11)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('Wye [m/s]')
    xlabel('time [s]')
    
%% plot wind velocity Wze
t = FilteredData{1}.State.Wze.time;
value = FilteredData{1}.State.Wze.data;
bounds = FilteredData{1}.State.Wze.std;

subplot(4,3,12)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('Wze [m/s]')
    xlabel('time [s]')      

%% Save figure
print('State1','-depsc')

%%
figure('Position', [40, 40, 800, 800]);

% Build plot title
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'State 2/3 (dark blue: state, grey dash: \pm1\sigma confidence bound)', ...
  'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', ...
  'Center', 'VerticalAlignment', 'Bottom' ) ;

%% plot Caup
t = FilteredData{1}.State.Caup.time;
value = FilteredData{1}.State.Caup.data;
bounds = FilteredData{1}.State.Caup.std;

subplot(5,2,1)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('Caup [rad]')
    xlabel('time [s]') 
     
%% plot Ca0
t = FilteredData{1}.State.Ca0.time;
value = FilteredData{1}.State.Ca0.data;
bounds = FilteredData{1}.State.Ca0.std;

subplot(5,2,2)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('Ca0 [-]')
    xlabel('time [s]')    
%% plot Cbsi
t = FilteredData{1}.State.Cbsi.time;
value = FilteredData{1}.State.Cbsi.data;
bounds = FilteredData{1}.State.Cbsi.std;

subplot(5,2,3)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('Cbsi [rad]')
    xlabel('time [s]')
   
%% plot Cb0
t = FilteredData{1}.State.Cb0.time;
value = FilteredData{1}.State.Cb0.data;
bounds = FilteredData{1}.State.Cb0.std;

subplot(5,2,4)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('Cb0 [-]')
    xlabel('time [s]')
    
%% plot lambda_Ax
t = FilteredData{1}.State.lambda_Ax.time;
value = FilteredData{1}.State.lambda_Ax.data;
bounds = FilteredData{1}.State.lambda_Ax.std;

subplot(5,2,5)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('\lambda_{Ax} [m/s^2]')
    xlabel('time [s]') 
    ylim([-0.2, 0.2])
    
%% plot lambda_Ay
t = FilteredData{1}.State.lambda_Ay.time;
value = FilteredData{1}.State.lambda_Ay.data;
bounds = FilteredData{1}.State.lambda_Ay.std;

subplot(5,2,6)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('\lambda_{Ay} [m/s^2]')
    xlabel('time [s]')    
    ylim([-0.2, 0.2])
    
%% plot lambda_Az
t = FilteredData{1}.State.lambda_Az.time;
value = FilteredData{1}.State.lambda_Az.data;
bounds = FilteredData{1}.State.lambda_Az.std;

subplot(5,2,7)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('\lambda_{Az} [m/s^2]')
    xlabel('time [s]')
    ylim([-0.2, 0.2])
    
%% plot lambda_p
t = FilteredData{1}.State.lambda_p.time;
value = FilteredData{1}.State.lambda_p.data;
bounds = FilteredData{1}.State.lambda_p.std;

subplot(5,2,8)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('\lambda_{p} [rad/s]')
    xlabel('time [s]')  
    ylim([-0.02, 0.02])
    
%% plot lambda_q
t = FilteredData{1}.State.lambda_q.time;
value = FilteredData{1}.State.lambda_q.data;
bounds = FilteredData{1}.State.lambda_q.std;

subplot(5,2,9)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('\lambda_{q} [rad/s]')
    xlabel('time [s]') 
    ylim([-0.02, 0.02])
    
%% plot lambda_r
t = FilteredData{1}.State.lambda_r.time;
value = FilteredData{1}.State.lambda_r.data;
bounds = FilteredData{1}.State.lambda_r.std;

subplot(5,2,10)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('\lambda_{r} [rad/s]')
    xlabel('time [s]')  
    ylim([-0.02, 0.02])
    
%% Save figure
print('State2','-depsc')



%%
figure('Position', [30, 30, 800, 800]);
sigma = 1;

% Build plot title
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'State 3/3 (dark blue: state, grey dash: \pm1\sigma confidence bound)', ...
  'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', ...
  'Center', 'VerticalAlignment', 'Bottom' ) ;

clr_est = [0.1 0.1 0.6];
clr_conf = [0.5 0.5 0.5]; 
         
%% plot body velocity u_L
t = FilteredData{1}.State.u_L.time;
value = FilteredData{1}.State.u_L.data;
bounds = FilteredData{1}.State.u_L.std;

subplot(2,3,1)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('u_L [m/s]')
    xlabel('time [s]')

%% plot body velocity u_R
t = FilteredData{1}.State.u_R.time;
value = FilteredData{1}.State.u_R.data;
bounds = FilteredData{1}.State.u_R.std;

subplot(2,3,2)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('u_R [m/s]')
    xlabel('time [s]')
    
%% plot body velocity v_L
t = FilteredData{1}.State.v_L.time;
value = FilteredData{1}.State.v_L.data;
bounds = FilteredData{1}.State.v_L.std;

subplot(2,3,3)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('v_L [m/s]')
    xlabel('time [s]')

%% plot body velocity v_R
t = FilteredData{1}.State.v_R.time;
value = FilteredData{1}.State.v_R.data;
bounds = FilteredData{1}.State.v_R.std;

subplot(2,3,4)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('v_R [m/s]')
    xlabel('time [s]')
    
%% plot body velocity w_L
t = FilteredData{1}.State.w_L.time;
value = FilteredData{1}.State.w_L.data;
bounds = FilteredData{1}.State.w_L.std;

subplot(2,3,5)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('w_L [m/s]')
    xlabel('time [s]')

%% plot body velocity w_R
t = FilteredData{1}.State.w_R.time;
value = FilteredData{1}.State.w_R.data;
bounds = FilteredData{1}.State.w_R.std;

subplot(2,3,6)
hold on
    grid on
    plot(t, (1 + sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, (1 - sigma*bounds).*value , '--', 'Color', clr_conf)
    plot(t, value, 'Color', clr_est)
    ylabel('w_R [m/s]')
    xlabel('time [s]')
    

%% Save figure
print('State3','-depsc')
end