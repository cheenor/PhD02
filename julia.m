function y=julia(theyear,themonth,theday); 
%return the julian day (year,day)定义函数julia为转�?�为julian日期函数，�?��?y为返回的julian日期，输入�?��?为theyear/年，themonth/月，theday/日 
leapyear=rem(theyear,4);    %年对4求余数 
leap=0;
if leapyear>0           
    leap=0;               %判断余数大于0，�?能被4整除，则用leap标记这一年�?是闰年（标记为0） 
else 
    if rem(theyear,100)==0 & rem(theyear,400)~=0 
        leap=0;           %能被4整除但是�?能被400整除也�?是闰年，则用leap标记这一年�?是闰年     
    else 
        leap=1;           %其他情况是闰年    
    end 
end
%%%%%%采用平年的日历，按1-12月分月计算julian日期    
   if themonth==1 
        juliaday=theday;     
   end 
    if themonth==2         
    juliaday=theday+31;     
   end 
    if themonth==3
    juliaday=theday+59;    
   end 
    if themonth==4 
        juliaday=theday+90;    
   end 
    if themonth==5 
        juliaday=theday+120;     
   end 
    if themonth==6 
        juliaday=theday+151;    
   end 
    if themonth==7 
        juliaday=theday+181;    
   end 
    if themonth==8 
        juliaday=theday+212;   
   end 
    if themonth==9 
        juliaday=theday+243;     
   end 
    if themonth==10         
       juliaday=theday+273;   
   end 
    if themonth==11 
        juliaday=theday+304;    
   end 
    if themonth==12       
       juliaday=theday+334;    
    end 
   if leap==1
       if themonth<=2            % 如果是闰年，当月份�?于2月时，julian日期与平年相�?�      
        juliaday=juliaday;    
       end 
       if themonth>=3 
        juliaday=juliaday+1;    %当月份大于2月时，在平年的基础上加一天    
       end 
    end 
y(1)=theyear;                 
y(2)=juliaday;