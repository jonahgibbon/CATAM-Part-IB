n=[10,30,50];
theta=2.2;
Counter=1;
dimN=size(n,2);

MLEStats=zeros(dimN,3);

while Counter<=dimN
    %Calculates random variables for different sample sizes
    data=zeros(n(Counter),200);
    mle_data=zeros(1,200);
    count=1;
    while count<=200
        u=rand(n(Counter),2);
        x=-log(1-u(:,1))/theta-log(1-u(:,2))/theta;
        data(:,count)=x;
        mle_data(1,count)=2*n(Counter)/sum(x);
        count=count+1;
    end

    %Calculating Average and Varience in data
    MLEStats(Counter,1)=sum(mle_data)/200;
    var_mle=0;
    weight_var_mle=0;
    count=1;
    while count<=200
        var_mle=var_mle+(mle_data(1,count)-MLEStats(Counter,1))^2;
        weight_var_mle=weight_var_mle+(mle_data(1,count)-theta)^2;
        count=count+1;
    end
    MLEStats(Counter,2)=var_mle/200;
    MLEStats(Counter,3)=weight_var_mle/200;
    
    if Counter==1
        XLower=min(mle_data)-0.2;
        XHigher=max(mle_data)+0.2;
    end
    
    %Plots histogram with appropiate data
    figure
    histogram(mle_data,'Normalization','pdf')
    xlim([XLower,XHigher])
    hold on 
    line([theta,theta],ylim,'LineWidth', 2,...
        'Color', 'r');
    hold on
    x=linspace(XLower,XHigher);
    y=pdf(x,2*n(Counter),theta);
    line(x,y,'LineWidth',2,'Color','g')
    
    legend('Data','\theta_{0}=2.2','Exact p.d.f','Location','northeast')
    xlabel('M.L.E')
    ylabel('Frequency Density')
    print(strcat('Image_8_',num2str(Counter)),'-depsc')
    
    Counter=Counter+1;
end
%MLE denotes average mle | Variance from average mle | Variance from 2.2
disp(MLEStats)
latex(sym(vpa(MLEStats)))

function answer = pdf(x,a,b)
    answer=(a*b)^(a)./((gamma(a)).*(x.^(a+1))).*exp(-b.*a./x);
end