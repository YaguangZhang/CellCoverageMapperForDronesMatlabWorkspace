lambda = 0.33; 
hr4=4;
hr3=3;
hr2=2;
ht=21;

axis=[];
p100=[];
p30=[];
p2=[];
pfsl=[];


for i=1000:5000
 d=10^(i/1000);
 axis =[axis d]; 
 fspower  = (lambda/(4*3.1415*d))^2 ;
 power100 = fspower * 4 *(sin(2*3.1415*ht*hr4/(lambda*d)))^2;
 power30  = fspower* 4 *(sin(2*3.1415*ht*hr3/(lambda*d)))^2;
 power2   = fspower * 4 *(sin(2*3.1415*ht*hr2/(lambda*d)))^2;

 p100 =[p100, 10*log10(power100)];
 p30 =[p30, 10*log10(power30)];
 p2 =[p2, 10*log10(power2)];
 pfsl=[pfsl, 10*log10(fspower)];
end

text('FontSize',18)

semilogx(axis,p100, 'g-',axis,p30, 'b-',axis,p2, 'r-',axis,pfsl,'y-')
 
xlabel('distance in m');
ylabel('pathloss');
text(1000,-66,'blue  : hr=30m');
text(1000,-74,'red   : hr=2m');
text(1000,-58,'red   : hr=100m');
text(1000,-50,'yellow: free space');

text(50,-180,'lambda = 0.33 m');
text(50,-190,'hr = 2 m');

grid on; grid minor;