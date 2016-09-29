base_filename = '080216-slm-pockels-x%s-y%s-r5-wffilter.txt'

y_range = -140:35:140;
x_range = y_range;
voltages = 0:.1:2;

power_curves = zeros(length(voltages),length(x_range),length(y_range));

figure

for x = 1:length(x_range)
    for y = 1:length(y_range)
        
        fullname = sprintf(base_filename,num2str(x_range(x)),num2str(y_range(y)));
        data = textread(fullname, '%f','headerlines',8);
        data = reshape(data,[2 length(voltages)])';
        power_curves(:,x,y) = data(:,2);
        
        subplot(length(y_range),length(x_range),(x-1)*length(y_range) + y)
        plot(data(:,1),data(:,2))
        
        ylim([0 250])
        
    end
end

%%
image1 = zeros(9,9);
image2 = zeros(9,9);

for x = 1:length(x_range)
    for y = 1:length(y_range)
        
        image1(x,y) = power_curves(end,x,y);
        image2(x,y) = power_curves(10,x,y);
        
    end
end

figure;
imagesc(image1)
figure; imagesc(image2)
figure; imagesc(image1./image2)

%%

voltage_measurments = 0:.1:2.0;
voltages = 0:.01:2.0;
num_volts = length(voltages);

power_curves_fit = zeros(num_volts,9,9);

for x = 1:length(x_range)
    for y = 1:length(y_range)
        
        pockels_fit = fit_pockels_curves(voltage_measurments,squeeze(power_curves(:,x,y))');
        power_curves_fit(:,x,y) = sin_sq(voltages,pockels_fit);
        
    end
end

%%
figure

for x = 1:length(x_range)
    for y = 1:length(y_range)
        
        
        subplot(length(y_range),length(x_range),(x-1)*length(y_range) + y)
        plot(voltages,power_curves_fit(:,x,y))
        
        ylim([0 250])
        
    end
end

%%

target_wattage = 10;
pockels_in = zeros(9,9);
power_predict = zeros(9,9);

for x = 1:length(x_range)
    for y = 1:length(y_range)
        
        target_ind = find(power_curves_fit(:,x,y) > target_wattage,1,'first');
        pockels_in(x,y) = voltages(target_ind);
        power_predict(x,y) = power_curves_fit(target_ind,x,y);
    end
end

figure;
imagesc(pockels_in)

figure;
imagesc(power_predict)







