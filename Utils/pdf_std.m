function stand_dev=pdf_std(prob_func, values) 
    y_approx=values.*prob_func; y_approx2=(values.^2).*prob_func;
    stand_dev=sqrt(sum(y_approx2)-sum(y_approx).^2);
end