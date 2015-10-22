figure

for j = 1:size(pockels_fits,1)
    
    i = j;
    ceil(i/7) + (mod(i-1,7))*7
    subplot(7,7,ceil(i/7) + (mod(i-1,7))*7)
    plot(x_volt,pockels_fits(i,1)*sin(pockels_fits(i,2)*x_volt + pockels_fits(i,3)).^2,'g',[50 100 125 150 180 199],roi_curve_measurements(:,i),'b+')
    ylim([0 280])
end