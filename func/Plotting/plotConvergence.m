function [] = plotConvergence(FilteredData)
%
figure('Position', [2500, -10, 210*3.8, 297*3.6]);
title('Convergence check for multiple initial conditions (1/2)')
Tmin = 0;
Tmax = 5;
Nlimit = 1;
%% plot x-position
t = FilteredData{1}.State.x.time;
value1 = FilteredData{1}.State.x.data;
value2 = FilteredData{2}.State.x.data;
value3 = FilteredData{3}.State.x.data;

subplot(4,3,1)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('x-position [m]')
    xlabel('time [s]') 
%% plot y-position
t = FilteredData{1}.State.y.time;
value1 = FilteredData{1}.State.y.data;
value2 = FilteredData{2}.State.y.data;
value3 = FilteredData{3}.State.y.data;

subplot(4,3,2)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('y-position [m]')
    xlabel('time [s]')   
%% plot z-position
t = FilteredData{1}.State.z.time;
value1 = FilteredData{1}.State.z.data;
value2 = FilteredData{2}.State.z.data;
value3 = FilteredData{3}.State.z.data;

subplot(4,3,3)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('z-position [m]')
    xlabel('time [s]')     
%% plot body velocity u
t = FilteredData{1}.State.u.time;
value1 = FilteredData{1}.State.u.data;
value2 = FilteredData{2}.State.u.data;
value3 = FilteredData{3}.State.u.data;

subplot(4,3,4)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('u [m/s]')
    xlabel('time [s]') 
%% plot body velocity v
t = FilteredData{1}.State.v.time;
value1 = FilteredData{1}.State.v.data;
value2 = FilteredData{2}.State.v.data;
value3 = FilteredData{3}.State.v.data;

subplot(4,3,5)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('v [m/s]')
    xlabel('time [s]')
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end   
%% plot body velocity w
t = FilteredData{1}.State.w.time;
value1 = FilteredData{1}.State.w.data;
value2 = FilteredData{2}.State.w.data;
value3 = FilteredData{3}.State.w.data;

subplot(4,3,6)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('w [m/s]')
    xlabel('time [s]')
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end 
%% plot euler angle phi
t = FilteredData{1}.State.phi.time;
value1 = FilteredData{1}.State.phi.data*180/pi;
value2 = FilteredData{2}.State.phi.data*180/pi;
value3 = FilteredData{3}.State.phi.data*180/pi;

subplot(4,3,7)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('phi [deg]')
    xlabel('time [s]')
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end 
%% plot euler angle theta
t = FilteredData{1}.State.theta.time;
value1 = FilteredData{1}.State.theta.data*180/pi;
value2 = FilteredData{2}.State.theta.data*180/pi;
value3 = FilteredData{3}.State.theta.data*180/pi;

subplot(4,3,8)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('theta [deg]')
    xlabel('time [s]')  
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end   
%% plot euler angle psi
t = FilteredData{1}.State.psi.time;
value1 = FilteredData{1}.State.psi.data*180/pi;
value2 = FilteredData{2}.State.psi.data*180/pi;
value3 = FilteredData{3}.State.psi.data*180/pi;

subplot(4,3,9)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('psi [deg]')
    xlabel('time [s]')       
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end
%% plot wind velocity Wxe
t = FilteredData{1}.State.Wxe.time;
value1 = FilteredData{1}.State.Wxe.data;
value2 = FilteredData{2}.State.Wxe.data;
value3 = FilteredData{3}.State.Wxe.data;

subplot(4,3,10)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('Wxe [m/s]')
    xlabel('time [s]')
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end 
%% plot wind velocity Wye
t = FilteredData{1}.State.Wye.time;
value1 = FilteredData{1}.State.Wye.data;
value2 = FilteredData{2}.State.Wye.data;
value3 = FilteredData{3}.State.Wye.data;

subplot(4,3,11)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('Wye [m/s]')
    xlabel('time [s]')
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end   
%% plot wind velocity Wze
t = FilteredData{1}.State.Wze.time;
value1 = FilteredData{1}.State.Wze.data;
value2 = FilteredData{2}.State.Wze.data;
value3 = FilteredData{3}.State.Wze.data;

subplot(4,3,12)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('Wze [m/s]')
    xlabel('time [s]')      
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end
%%
figure('Position', [2500, -10, 210*3.8, 297*3.6]);
title('Convergence check for multiple initial conditions (2/2)')
%% plot Caup
t = FilteredData{1}.State.Caup.time;
value1 = FilteredData{1}.State.Caup.data;
value2 = FilteredData{2}.State.Caup.data;
value3 = FilteredData{3}.State.Caup.data;

subplot(5,2,1)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('Caup [rad]')
    xlabel('time [s]') 
    
if isempty(t) ~= 1    
    checkConverge = abs((abs(value1 - value2) + abs(value1 - value3))/mean(value1));
    Nconverge = checkConverge < Nlimit; idx = find(Nconverge~=0, 1, 'first');
    Tmin = 0; Tmax = t(idx);
    xlabel('time [s]')
    xlim([Tmin Tmax])
end       
%% plot Ca0
t = FilteredData{1}.State.Ca0.time;
value1 = FilteredData{1}.State.Ca0.data;
value2 = FilteredData{2}.State.Ca0.data;
value3 = FilteredData{3}.State.Ca0.data;

subplot(5,2,2)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('Ca0 [-]')
    xlabel('time [s]')
    
if isempty(t) ~= 1    
    checkConverge = abs((abs(value1 - value2) + abs(value1 - value3))/mean(value1));
    Nconverge = checkConverge < Nlimit; idx = find(Nconverge~=0, 1, 'first');
    Tmin = 0; Tmax = t(idx);
    xlabel('time [s]')
    xlim([Tmin Tmax])
end           
%% plot Cbsi
t = FilteredData{1}.State.Cbsi.time;
value1 = FilteredData{1}.State.Cbsi.data;
value2 = FilteredData{2}.State.Cbsi.data;
value3 = FilteredData{3}.State.Cbsi.data;

subplot(5,2,3)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('Cbsi [rad]')
    xlabel('time [s]')
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end
%% plot Cb0
t = FilteredData{1}.State.Cb0.time;
value1 = FilteredData{1}.State.Cb0.data;
value2 = FilteredData{2}.State.Cb0.data;
value3 = FilteredData{3}.State.Cb0.data;

subplot(5,2,4)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('Cb0 [-]')
    xlabel('time [s]')
    if isempty(t) ~= 1   
        xlim([Tmin Tmax])    
    end   
%% plot lambda_Ax
t = FilteredData{1}.State.lambda_Ax.time;
value1 = FilteredData{1}.State.lambda_Ax.data;
value2 = FilteredData{2}.State.lambda_Ax.data;
value3 = FilteredData{3}.State.lambda_Ax.data;
subplot(5,2,5)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('\lambda_{Ax} [m/s^2]')
    xlabel('time [s]') 
if isempty(t) ~= 1    
    checkConverge = abs((abs(value1 - value2) + abs(value1 - value3))/mean(value1));
    Nconverge = checkConverge < Nlimit; idx = find(Nconverge~=0, 1, 'first');
    Tmin = 0; Tmax = t(idx);
    xlabel('time [s]')
    xlim([Tmin Tmax])
end         
%% plot lambda_Ay
t = FilteredData{1}.State.lambda_Ay.time;
value1 = FilteredData{1}.State.lambda_Ay.data;
value2 = FilteredData{2}.State.lambda_Ay.data;
value3 = FilteredData{3}.State.lambda_Ay.data;

subplot(5,2,6)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('\lambda_{Ay} [m/s^2]')
    xlabel('time [s]')    
if isempty(t) ~= 1    
    checkConverge = abs((abs(value1 - value2) + abs(value1 - value3))/mean(value1));
    Nconverge = checkConverge < Nlimit; idx = find(Nconverge~=0, 1, 'first');
    Tmin = 0; Tmax = t(idx);
    xlabel('time [s]')
    xlim([Tmin Tmax])
end       
%% plot lambda_Az
t = FilteredData{1}.State.lambda_Az.time;
value1 = FilteredData{1}.State.lambda_Az.data;
value2 = FilteredData{2}.State.lambda_Az.data;
value3 = FilteredData{3}.State.lambda_Az.data;

subplot(5,2,7)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('\lambda_{Az} [m/s^2]')
    xlabel('time [s]')
if isempty(t) ~= 1    
    checkConverge = abs((abs(value1 - value2) + abs(value1 - value3))/mean(value1));
    Nconverge = checkConverge < Nlimit; idx = find(Nconverge~=0, 1, 'first');
    Tmin = 0; Tmax = t(idx);
    xlabel('time [s]')
    xlim([Tmin Tmax])
end          
%% plot lambda_p
t = FilteredData{1}.State.lambda_p.time;
value1 = FilteredData{1}.State.lambda_p.data;
value2 = FilteredData{2}.State.lambda_p.data;
value3 = FilteredData{3}.State.lambda_p.data;

subplot(5,2,8)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('\lambda_{p} [rad/s]')
    xlabel('time [s]')  
if isempty(t) ~= 1    
    checkConverge = abs((abs(value1 - value2) + abs(value1 - value3))/mean(value1));
    Nconverge = checkConverge < Nlimit; idx = find(Nconverge~=0, 1, 'first');
    Tmin = 0; Tmax = t(idx);
    xlabel('time [s]')
    xlim([Tmin Tmax])
end       
%% plot lambda_q
t = FilteredData{1}.State.lambda_q.time;
value1 = FilteredData{1}.State.lambda_q.data;
value2 = FilteredData{2}.State.lambda_q.data;
value3 = FilteredData{3}.State.lambda_q.data;

subplot(5,2,9)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('\lambda_{q} [rad/s]')
    xlabel('time [s]') 
if isempty(t) ~= 1    
    checkConverge = abs((abs(value1 - value2) + abs(value1 - value3))/mean(value1));
    Nconverge = checkConverge < Nlimit; idx = find(Nconverge~=0, 1, 'first');
    Tmin = 0; Tmax = t(idx);
    xlabel('time [s]')
    xlim([Tmin Tmax])
end     
%% plot lambda_r
t = FilteredData{1}.State.lambda_r.time;
value1 = FilteredData{1}.State.lambda_r.data;
value2 = FilteredData{2}.State.lambda_r.data;
value3 = FilteredData{3}.State.lambda_r.data;

subplot(5,2,10)
hold on
    plot(t, value1, 'r')
    plot(t, value2 ,'g--')
    plot(t, value3 ,'g--')
    ylabel('\lambda_{r} [rad/s]')
    xlabel('time [s]')   
if isempty(t) ~= 1    
    checkConverge = abs((abs(value1 - value2) + abs(value1 - value3))/mean(value1));
    Nconverge = checkConverge < Nlimit; idx = find(Nconverge~=0, 1, 'first');
    Tmin = 0; Tmax = t(idx);
    xlabel('time [s]')
    xlim([Tmin Tmax])
end      

%% Save figure
print('Convergence','-depsc')

end