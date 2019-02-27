close all
clear all

dt = 0.0007717 ;
time = dt:dt:(dt*1000) ; 
path = [ '~/Documents/lab/KD-project/AdvectionDiffusion/right/asi2/p2/RT1/' ] ;

    
volume_integral = load( [path , 'int_r.dat'] );
volume = load( [path,'Volume.dat'] ) ;

figure
plot(volume_integral(:,end)) ;

l=length(volume_integral(:,1)) ;

time_integral = trapz(time ,volume_integral(l-1000 + 1:end,2) ) ;
RT1 = time_integral/volume ;

norm_volume_integral = volume_integral/volume - 1 ;
norm_RT1 = trapz(time ,norm_volume_integral(l-1000 + 1:end,2) ) ;

display(path)
display(['RT1 = ', num2str(norm_RT1), ' [s]'])
display(['Vol = ', num2str(volume),])