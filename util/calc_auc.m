function auc = calc_auc(x, do_plot, label, show_title)
errors = x;
errors = sort(errors);
errors(errors > .25) = .25;
prec = [errors; 0.25];
rec = [1:length(errors) length(errors)]./(length(errors));
auc = trapz(prec, rec)/0.25;
if do_plot
plot(prec, rec, '.-')
if show_title
    title(strcat(label, sprintf(': %2.4f', auc*100)))
end
disp(sprintf('AUC: %2.4f', auc*100))
grid on
else
    fprintf(strcat(label, sprintf(': %2.4f', auc*100),'\n'));
end
