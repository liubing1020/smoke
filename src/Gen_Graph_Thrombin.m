clear variables
addpath ode/
addpath util/
addpath analysis/
addpath ../models/thrombin/

cd ode
compile('thrombin','graph')
cd ..
model = thrombin_pathway();

%(timse discretization)
dt = 50;
% Number of discrete time steps
nTime = 20;

[numData,textData,rawData]=xlsread('../models/thrombin/training.csv');
[numData2,textData2,rawData2]=xlsread('../models/thrombin/test.csv');

%p=[6.0865200E+00	2.6745810E+01	4.6760910E+00	1.6768470E-02	1.3037070E-03	9.2704960E-01	6.8456630E-03	6.0252380E-02	1.0650330E-02	4.1131490E+01	2.4171830E+01	1.9003880E-01	1.1071550E-02	6.2257870E-02	2.3942280E+00	9.9769990E-01	2.2243080E+01	4.9643860E-02	1.0056920E-02	1.1666960E-02	6.1986880E-04	1.9507600E+01	6.6597860E-01	1.7504320E+01	1.1554470E+00	1.1396470E-02	1.7183000E-02	2.3423680E-03	3.9566250E+01	1.8928750E+00	3.8276150E+00	1.0451280E+00	3.1291330E+01	7.9217210E-01	1.8091650E+00	1.8058840E+00	1.4974130E-02	1.5214510E+01	5.5516200E+00	7.6531780E+00	8.2715960E+00	6.9473880E+01	2.7018580E-01	3.8760780E+00	1.6358850E-01	6.1092330E+00	1.0068980E+02	9.2853450E+00	9.3989770E+01	7.9304970E+00	6.7436380E+01	1.2540300E+01	2.7898170E+02	1.0813390E+01	2.1465810E+02	2.5945570E-04	6.9487440E-02	5.2467750E-01	8.8996010E-01	5.8999310E+01	5.1773430E+00	7.6859950E-02	7.6937960E+00	3.0876680E+00	4.7649290E+00	7.3994230E+00	5.5929790E+00	4.7813580E-01	1.3000580E+02	8.8673350E-01	7.4697630E+00	2.1257170E+00	1.6583370E+01	1.1346980E+00	9.2773250E-01	1.4752810E+01	1.7709650E+00	8.7929610E-02	2.4504790E+01	5.4547510E-02	2.8728120E+01	8.4802100E+01	2.7370510E+00	1.7825430E+01	1.4399250E+01	1.3274500E+01	5.4149930E+00	6.7851310E+01	3.4647580E+00	5.2176780E+01	3.2393180E+00	1.6298120E+00	3.9087260E-01	6.4584300E+00	5.1986880E+00	5.0378950E+00	3.3432820E+00	4.3034080E+00	8.7603050E-03	1.6152200E-02];
%p=[7.086273e+00	2.846783e+01	2.993835e+00	1.002802e-02	3.757101e-03	1.028638e+00	7.717726e-03	9.798432e-02	6.015651e-03	3.958323e+01	3.572978e+01	1.058208e-01	8.696442e-03	4.860755e-02	2.371172e+00	1.453982e+00	1.397023e+01	4.979861e-02	1.373801e-02	1.160183e-02	5.859644e-04	1.399058e+01	7.450388e-01	2.133391e+01	9.679932e-01	5.421983e-03	1.274266e-02	3.861512e-03	2.910297e+01	1.911989e+00	1.475083e+00	6.310438e-01	4.218668e+01	1.882034e+00	3.728079e+00	1.546253e+00	8.360688e-03	2.127438e+01	1.428831e+01	2.862680e+00	1.182315e+01	4.844146e+01	2.863678e-01	2.601233e+00	1.910602e-01	1.934243e+01	2.517347e+02	1.271577e+01	6.288036e+01	1.350696e+01	6.094905e+01	1.064768e+01	5.187407e+02	1.957556e+01	2.660721e+02	3.264372e-04	1.462063e-01	3.144535e-01	9.522856e-01	3.777836e+01	4.500517e+00	7.336448e-02	5.190382e+00	6.526317e+00	6.911271e+00	4.508638e+00	4.237464e+00	6.462450e-01	1.112125e+02	8.613193e-01	5.608215e+00	2.297309e+00	1.417417e+01	1.563710e+00	2.096827e+00	1.503588e+01	2.242330e+00	1.333588e-01	3.316962e+01	1.607059e-02	3.025983e+01	8.966221e+01	3.201297e+00	1.565761e+01	2.153539e+01	1.373546e+01	9.945317e+00	1.069430e+02	2.399135e+00	4.400815e+01	3.536527e+00	1.824386e+00	3.254307e-01	3.739225e+00	3.874536e+00	2.881284e+00	5.692452e+00	6.691374e+00	7.779338e-03	1.404319e-02];
%0.8-p
%p=[5.79567E+00	1.37325E+01	2.24985E+00	8.29185E-03	1.29457E-03	6.92517E-01	7.31094E-03	4.69415E-02	1.83317E-02	4.53349E+01	1.33487E+01	1.82316E-01	1.10912E-02	2.63494E-02	3.78661E+00	1.15677E+00	1.07212E+01	3.18870E-02	6.41022E-03	1.16965E-02	5.19047E-04	5.98109E+00	5.33486E-01	2.92919E+01	1.27356E+00	7.13426E-03	2.45780E-02	1.96704E-03	3.49066E+01	1.68656E+00	1.37353E+00	1.77805E+00	5.71653E+01	1.26128E+00	5.48163E+00	1.72638E+00	9.72328E-03	3.90971E+01	1.69140E+01	3.62757E+00	1.97851E+01	9.25645E+01	1.73213E-01	3.80057E+00	1.83822E-01	1.40752E+01	1.02654E+02	1.26996E+01	8.96704E+01	9.20554E+00	5.32755E+01	9.20258E+00	5.23694E+02	1.01679E+01	2.60170E+02	4.94275E-04	1.51285E-01	2.44448E-01	7.01808E-01	2.18953E+01	2.29315E+00	7.72713E-02	5.74857E+00	2.19929E+00	3.57702E+00	6.59939E+00	3.15795E+00	7.10237E-01	1.31088E+02	8.01110E-01	8.02679E+00	1.81503E+00	2.76327E+01	1.81190E+00	1.70948E+00	1.55108E+01	2.74358E+00	9.33356E-02	3.46228E+01	5.07006E-02	3.40492E+01	1.13148E+02	2.59533E+00	2.62188E+01	1.47070E+01	1.23903E+01	9.02117E+00	1.11937E+02	1.77787E+00	3.98241E+01	3.62119E+00	3.43844E+00	2.93813E-01	4.66430E+00	6.27566E+00	4.93716E+00	5.49062E+00	5.51718E+00	7.88569E-03	8.69123E-03];
%0.8
%p=[6.82566E+00	9.54922E+00	6.33851E+00	8.59988E-03	2.34000E-03	8.98617E-01	1.07602E-02	5.29924E-02	1.95238E-02	2.13872E+01	3.54643E+01	1.61884E-01	2.03193E-02	9.00257E-02	2.07253E+00	1.40784E+00	1.58502E+01	6.52682E-02	7.99855E-03	1.17784E-02	4.22263E-04	1.65282E+01	9.28609E-01	2.61718E+01	1.07490E+00	6.27419E-03	2.21698E-02	4.04330E-03	3.96437E+01	1.39078E+00	3.27869E+00	1.38773E+00	3.37204E+01	1.05281E+00	4.24665E+00	1.47158E+00	1.15715E-02	1.58436E+01	5.94203E+00	3.48346E+00	8.00455E+00	7.02765E+01	2.38629E-01	3.50424E+00	1.27343E-01	8.26929E+00	1.84797E+02	1.05728E+01	6.45777E+01	6.88439E+00	2.68711E+01	1.27202E+01	6.84019E+02	9.96194E+00	1.21795E+02	4.73605E-04	1.75656E-01	2.44737E-01	7.77791E-01	2.19327E+01	4.45766E+00	6.21505E-02	5.09580E+00	7.38866E+00	7.64410E+00	6.65278E+00	6.85209E+00	3.11532E-01	8.72749E+01	4.79324E-01	6.91449E+00	6.51479E-01	2.49446E+01	2.36212E+00	1.95238E+00	1.29908E+01	2.28064E+00	9.59235E-02	3.36345E+01	4.92752E-02	2.96618E+01	1.06945E+02	2.58721E+00	2.16760E+01	1.13470E+01	1.30336E+01	1.15943E+01	4.37203E+01	2.43794E+00	4.26891E+01	2.69224E+00	1.84530E+00	3.67310E-01	4.73517E+00	4.71843E+00	4.84140E+00	4.66242E+00	4.49967E+00	8.39030E-03	8.94971E-03];
%0.9
p=[7.70122E+00	2.63813E+01	2.49982E+00	1.29956E-02	2.30501E-03	7.41003E-01	3.68614E-03	4.55653E-02	1.28817E-02	3.36550E+01	1.13878E+01	1.03366E-01	1.12431E-02	7.84641E-02	2.03662E+00	1.27471E+00	1.76889E+01	5.64603E-02	1.07996E-02	1.18202E-02	5.49217E-04	1.51114E+01	5.03952E-01	2.42941E+01	1.21310E+00	6.49747E-03	2.44241E-02	2.15543E-03	4.20206E+01	9.48793E-01	3.35230E+00	7.51696E-01	5.54628E+01	1.53948E+00	2.17517E+00	1.79012E+00	2.28697E-02	1.46871E+01	1.79529E+01	5.74335E+00	1.29397E+01	7.18005E+01	1.73436E-01	2.95089E+00	1.55444E-01	1.86071E+01	1.02042E+02	8.57435E+00	6.75978E+01	1.09272E+01	4.65044E+01	1.51767E+01	5.56979E+02	6.76375E+00	1.71869E+02	5.10325E-04	1.74350E-01	3.34644E-01	6.46509E-01	3.16072E+01	1.92856E+00	5.82235E-02	3.94996E+00	6.73280E+00	3.70777E+00	6.95748E+00	4.17618E+00	4.35925E-01	1.22024E+02	8.30819E-01	2.33185E+00	1.63691E+00	9.26820E+00	2.33314E+00	1.96193E+00	1.42312E+01	2.29638E+00	1.62699E-01	2.96219E+01	3.35651E-02	2.65977E+01	3.79108E+01	2.64654E+00	3.60514E+01	1.38730E+01	9.56272E+00	9.39370E+00	5.94901E+01	1.57000E+00	4.16739E+01	1.41518E+00	2.07260E+00	2.81455E-01	6.93564E+00	3.76962E+00	4.86609E+00	4.85880E+00	3.74033E+00	9.26463E-03	7.21367E-03];
%0.9-p
%p=[4.85457E+00	2.87110E+01	4.14011E+00	1.37244E-02	1.28063E-03	8.05434E-01	5.96058E-03	8.86733E-02	1.99271E-02	2.79154E+01	2.92629E+01	1.56648E-01	1.31031E-02	6.86958E-02	2.20129E+00	1.08860E+00	9.87707E+00	6.69914E-02	8.22254E-03	1.17231E-02	6.49262E-04	8.65828E+00	3.39878E-01	1.58776E+01	1.12144E+00	5.11879E-03	1.04397E-02	3.23036E-03	1.48830E+01	9.01417E-01	1.72719E+00	1.93183E+00	3.71926E+01	6.62201E-01	1.84593E+00	1.60437E+00	1.22907E-02	3.21825E+01	1.34212E+01	2.82068E+00	1.16575E+01	7.44479E+01	2.08152E-01	3.82155E+00	1.89216E-01	7.25958E+00	1.84002E+02	1.54310E+01	5.29824E+01	8.06944E+00	7.41503E+01	5.61548E+00	6.71334E+02	1.07653E+01	1.49573E+02	5.26599E-04	1.16877E-01	3.90280E-01	7.85624E-01	3.24193E+01	5.10617E+00	5.28783E-02	4.71036E+00	5.65355E+00	2.37464E+00	7.01484E+00	4.85898E+00	8.39254E-01	1.07196E+02	5.39828E-01	4.94254E+00	1.86705E+00	1.57081E+01	8.97346E-01	1.99040E+00	8.45851E+00	1.30151E+00	1.68991E-01	2.67044E+01	2.06408E-02	3.20516E+01	1.10541E+02	3.77827E+00	4.29255E+01	1.05465E+01	1.68727E+01	9.56870E+00	4.19199E+01	1.72629E+00	5.25299E+01	1.08566E+00	1.93731E+00	2.11003E-01	6.65706E+00	5.61656E+00	4.49782E+00	4.89993E+00	2.93997E+00	9.29858E-03	1.78674E-02];




B = model.p_nominal;

for i=1:length(model.p_estimate)
	B(:,(model.p_estimate(i)))= p(:,i);
end

samples=1000;
for pp=1:samples
	%pRand=model.p_nominal;
	disp(pp);
	B = samplingScheme(B,'uniform',0.0025);
     for i=1:length(model.p_nominal)
             if find(model.p_estimate==i)
             else 
                 B(i)= model.p_nominal(i);
             end
     end
	xInit = samplingScheme(model.x_init,'uniform',0.025);
	tf = Gen_Graph(B,xInit,nTime,dt);
	AB{pp} = transpose(reshape(tf,length(model.x_init),nTime));
end

t= linspace(0,nTime*dt,nTime);

a=[6,16];
		
		
save('thrombin_paramest')	
%%
load('thrombin_paramest')
figure;
idx = [11 23 41 53 57 58 66 104];
chooseidx = [2 3 8];
xidx = idx(chooseidx);
nplot = length(xidx)+1;


title('Training');
for i=1:nplot-1
	subplot(1,nplot,i);		
	for pp=1:samples
		axis([0 nTime*dt 0 1.2*max(max(AB{pp}(:,xidx(i))),max(numData(:,chooseidx(i)+1)))]);
		plot(t, AB{pp}(:,xidx(i)),'Color',[0.5 0.5 0.5]);
		hold on;
	end
	
	title(strcat(model.x_names{xidx(i)}),'FontSize',13);

	errorbar((numData(:,1)*dt),numData(:,chooseidx(i)+1),numData(:,chooseidx(i)+1)*0.1,'ks','MarkerEdgeColor','k',...
	 'MarkerFaceColor','w','MarkerSize',5);
	 hold on;
	set(gca,'FontSize',12)
end

subplot(1,nplot,nplot);	
specid = 16;
for pp=1:samples
	plot(t,AB{pp}(:,specid),'Color',[0.5 0.5 0.5]);
	hold on;
end
set(gca,'FontSize',12)
title(strcat(model.x_names{specid}),'FontSize',13);

%%
figure;
specid = 86;
for pp=1:samples
	plot(t,AB{pp}(:,specid),'Color',[0.5 0.5 0.5]);
	hold on;
end
set(gca,'FontSize',12)
title(strcat(model.x_names{specid}),'FontSize',13);

errorbar((numData2(:,1)*dt),numData2(:,3),numData2(:,3)*0.1,'ks','MarkerEdgeColor','k',...
		 'MarkerFaceColor','w','MarkerSize',5);