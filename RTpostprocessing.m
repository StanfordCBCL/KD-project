close all
clear all

dt =  0.000811;
time = dt:dt:(dt*1000) ; 
path = [ '~/Documents/lab/KD-project/AdvectionDiffusion/right/asi2/d5/RT1/' ] ;
%path = ['~/Documents/lab/KD-project/AdvectionDiffusion/baseline/asi6/' ];

no_offset = true;

    
volume_integral = load( [path , 'int_r.dat'] );
volume = load( [path,'Volume.dat'] ) ;

figure
plot(volume_integral(:,end)) ;

l=length(volume_integral(:,1)) ;

% initialize variables
norm_RT1 = 0;
RT1 = 0;
volume = 0;

if no_offset
    time_integral = trapz(time ,volume_integral(l-1000 + 1:end,2) ) ;
    RT1 = time_integral/volume ;

    norm_volume_integral = volume_integral/volume - 1 ;
    norm_RT1 = trapz(time ,norm_volume_integral(l-1000 + 1:end,2) ) ;
    RT1 = trapz(time ,volume_integral(l-1000 + 1:end,2) )/volume 
    
else
    
    time_integral = trapz(time ,volume_integral(l-1000 + 1:end,2) ) ;
    RT1 = time_integral/volume ;

    norm_volume_integral = volume_integral/volume - 1 ;
    norm_RT1 = trapz(time ,norm_volume_integral(l-1000 + 1:end,2) ) ;

end

display(path)
display(['RT1 = ', num2str(norm_RT1), ' [s]'])
display(['Vol = ', num2str(volume),])
