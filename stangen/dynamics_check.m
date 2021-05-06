% #####################################################
%
% Linkage system: Linkage system of Ghassaei
% 
% Seppe Vilain, Arnoud Deboeck
% 
% #####################################################

function dynamics_check(F1_2_x,F2_12_x,F2_3_x,F3_4_x,F3_5_x,F4_6_x,F6_7_x,F6_8_x,F8_9_x, ...
          F8_10_x,F12_10_x,F12_11_x,F1_5_x,F1_7_x,F1_9_x,F1_11_x,t,plot_dynCheck)

        if plot_dynCheck
        figure

        subplot(8,1,1)
        plot(t,F1_2_x)
        xlabel('t [s]')
        ylabel('F1_2_x [m/s] ')
        axis tight

        subplot(8,1,2)
        plot(t,F2_3_x)
        xlabel('t [s]')
        ylabel('F2_3_x [m/s] ')
        axis tight

        subplot(8,1,3)
        plot(t,F2_12_x)
        xlabel('t [s]')
        ylabel('F2_12_x [m/s] ')
        axis tight

        subplot(8,1,4)
        plot(t,F3_4_x)
        xlabel('t [s]')
        ylabel('F3_4_x [m/s] ')
        axis tight

        subplot(8,1,5)
        plot(t,F3_5_x)
        xlabel('t [s]')
        ylabel('F1_2_x [m/s] ')
        axis tight

        subplot(8,1,6)
        plot(t,F4_6_x)
        xlabel('t [s]')
        ylabel('F4_6_x [m/s] ')
        axis tight

        subplot(8,1,7)
        plot(t,F6_7_x)
        xlabel('t [s]')
        ylabel('F6_7_x [m/s] ')
        axis tight

        subplot(8,1,8)
        plot(t,F6_8_x)
        xlabel('t [s]')
        ylabel('F6_8_x [m/s] ')
        axis tight
        end

        if plot_dynCheck
        figure

        subplot(8,1,1)
        plot(t,F8_9_x)
        xlabel('t [s]')
        ylabel('F8_9_x [m/s] ')
        axis tight

        subplot(8,1,2)
        plot(t,F8_10_x)
        xlabel('t [s]')
        ylabel('F8_10_x [m/s] ')
        axis tight

        subplot(8,1,3)
        plot(t,F12_10_x)
        xlabel('t [s]')
        ylabel('F12_10_x [m/s] ')
        axis tight

        subplot(8,1,4)
        plot(t,F12_11_x)
        xlabel('t [s]')
        ylabel('F12_11_x [m/s] ')
        axis tight

        subplot(8,1,5)
        plot(t,F1_5_x)
        xlabel('t [s]')
        ylabel('F1_5_x [m/s] ')
        axis tight

        subplot(8,1,6)
        plot(t,F1_7_x)
        xlabel('t [s]')
        ylabel('F1_7_x [m/s] ')
        axis tight

        subplot(8,1,7)
        plot(t,F1_9_x)
        xlabel('t [s]')
        ylabel('F1_9_x [m/s] ')
        axis tight

        subplot(8,1,8)
        plot(t,F1_11_x)
        xlabel('t [s]')
        ylabel('F1_2_x [m/s] ')
        axis tight
        end
end

