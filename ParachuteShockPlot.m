function [] = ParachuteShockPlot(dData, mData, tag, n, delay, infdelay, SItag)

% Plot a bunch of pretty important graphs from data passed from ParachuteShockCalculationFunction.m

if (contains(tag,"drogue",'IgnoreCase',true) == 1) || (contains(tag,"both",'IgnoreCase',true) == 1)

    figure(1)
    hold on; grid on;
    title('Drogue Parachute Shock vs. Time')
    plot(dData(:,1),dData(:,7));
    xlabel("Time (s)"); ylabel("Shock Force (Gs)")
    yyaxis right
    plot(dData(:,1),dData(:,6));
    xlabel("Time (s)"); 
    if (SItag == true)
        xlabel("Time (s)"); ylabel("Shock Force (N)");        
    else 
        xlabel("Time (s)"); ylabel("Shock Force (lbf)");
    end
    xline(delay,'--'); xline(delay + infdelay,'--')
    
    figure(2)
    title('Rocket Velocity at Apogee Deployment vs. Time')
    hold on; grid on;
    plot(dData(:,1),dData(:,3));
    xlabel("Time (s)"); 
    if (SItag == true)
        ylabel("Vertical Velocity (m/s)");
    else
        ylabel("Vertical Velocity (ft/s)");
    end
    xline(delay,'--'); xline(delay + infdelay,'--')

end

if (contains(tag,"main",'IgnoreCase',true) == 1) || (contains(tag,"both",'IgnoreCase',true) == 1)

    figure(n+1)
    title('Main Parachute Shock vs. Time')
    hold on; grid on;
    plot(mData(:,1),mData(:,7));
    xlabel("Time (s)"); ylabel("Shock Force (Gs)")
    yyaxis right
    plot(mData(:,1),mData(:,6));
    xlabel("Time (s)"); 
    if SItag == true
        ylabel("Shock Force (N)");     
    else 
        ylabel("Shock Force (lbf)");
    end
    xline(infdelay,'--');
    
    figure(n+2)
    title('Rocket Velocity at Main deployment vs. Time')
    hold on; grid on;
    plot(mData(:,1),mData(:,3));
    xlabel("Time (s)"); 
    if SItag == true
        ylabel("Vertical Velocity (m/s)");
    else 
        ylabel("Vertical Velocity (ft/s)");
    end
    xline(infdelay,'--');

end


end
