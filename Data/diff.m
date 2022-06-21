clc
close all
clear all

result = 'Result_';
true = 'True_';
list = {'S01_T01','S02_T01','S02_T02','S03_T02','S04_T02','S05_T02','S06_T01','S06_T02','S07_T02','S08_T01'};

for l = 1:size(list,2)
    d_result = load(strcat(result, list{l}));
    d_true = load(strcat(true, list{l}));
    d_result = d_result.BPM';
    d_true = d_true.BPM0;
    size(d_result);
    size(d_true);
    
    h = figure;
    plot(d_result);
    hold on;
    plot(d_true,'-r');
    hold off;
    title(list{l});
    xlabel('Time (sec)');
    ylabel('Heart Rate (BPM)');
    legend('prediction','true')
    av(l) = sum(abs(d_result-d_true))/size(d_result,1);
    disp(strcat(list{l},':', num2str(av(l))));
    saveas(h,strcat(list{l},'.jpg'));
end
av
sum(av)/length(av)