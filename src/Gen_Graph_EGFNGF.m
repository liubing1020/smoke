clear variables
addpath ode
addpath util
addpath analysis
addpath ../models/egf-ngf/

cd ode
compile('egfngf','graph')
cd ..
model = egf_ngf_pathway();

[numData,textData,rawData]=xlsread('../models/egf-ngf/training.csv');
[numData2,textData2,rawData2]=xlsread('../models/egf-ngf/test.csv');

        
% Time steps (timse discretization)
dt = 1;
% Number of discrete time steps
nTime = 61;

%p=[2.23207E-05	1.64508E-02	1.36877E-07	1.06305E-02	5.43478E+01	2.02448E+04	6.97616E-01	1.92567E+02	1.06946E+01	3.55225E-02	1.46373E+06	2.11166E+01	6.72084E-02	8.54658E+05	1.84173E+02	8.40492E+02	8.97220E-01	2.25511E+01	2.58758E+00	1.76656E+06];
%ga-pval
%p=[2.2320720E-05	1.6450780E-02	1.3687680E-07	1.0630470E-02	5.4347800E+01	2.0244820E+04	6.9761620E-01	1.9256700E+02	1.0694560E+01	3.5522480E-02	1.4637250E+06	2.1116550E+01	6.7208440E-02	8.5465790E+05	1.8417340E+02	8.4049200E+02	8.9721990E-01	2.2551050E+01	2.5875830E+00	1.7665580E+06];
%sres-pval
%p=[3.2295440E-05	8.1409570E-03	1.2778850E-07	6.3113630E-03	6.0548350E+01	4.1841460E+04	1.3347490E+00	7.6141900E+01	5.3480240E+00	3.5211490E-02	1.3211060E+06	1.7135650E+01	7.0503250E-02	8.7128140E+05	2.3299920E+02	8.4779940E+03	2.3159840E+00	3.5688490E+01	1.9835430E+00	8.8904490E+05];
%sres-normal
%p=[2.7302730E-05	1.1656730E-02	1.3380500E-07	7.8604320E-03	5.5419140E+01	3.7052400E+04	1.7173330E+00	4.7872900E+01	6.1632750E+00	3.5448820E-02	1.3479770E+06	1.3429600E+01	9.8553790E-02	1.2685280E+06	2.6680090E+02	1.3628370E+04	2.1307620E+00	3.4376380E+01	2.3683580E+00	1.2230730E+06];
%p=[2.352556e-05	1.029213e-02	1.368768e-07	1.063047e-02	5.600625e+01	1.934946e+04	6.976162e-01	2.845608e+02	1.089768e+01	2.335322e-02	9.933571e+05	2.099155e+01	6.720844e-02	8.606379e+05	2.336537e+02	1.901106e+03	8.347199e-01	2.305105e+01	3.587583e+00	1.991671e+06];
%sres 0.8-p
%p=[8.96327E-05	6.26203E-02	1.36640E-07	7.95575E-03	1.65925E+02	3.07178E+05	1.29793E+00	5.94778E+01	6.88517E+00	1.57824E-01	6.30486E+06	2.71575E+01	4.46288E-01	6.20015E+06	4.44714E+02	1.58366E+04	9.78970E+00	2.39416E+02	1.45775E+01	5.99365E+06];
%sres 0.8-n
%p=[9.69097E-05	1.15551E-02	1.35272E-07	8.14749E-03	4.91386E+01	3.27527E+05	2.20163E+00	7.70169E+01	1.36200E+01	1.62189E-01	6.28327E+06	1.22193E+01	4.35951E-01	5.86584E+06	3.85315E+02	2.82872E+04	1.85797E+00	4.00265E+01	4.90565E+00	3.74434E+06];
%sres 0.9p
%p=[1.47356E-04	9.10892E-03	1.37237E-07	9.02047E-03	2.53026E+02	1.11318E+05	8.48359E-01	4.74428E+01	5.03561E+00	1.92826E-01	7.60915E+06	8.66101E+01	4.35358E-01	6.23786E+06	8.06170E+02	2.70794E+04	1.03262E+01	2.19673E+02	1.89913E+01	6.77575E+06];

%sres 0.9n
p=[7.50669E-05	5.83982E-02	1.38100E-07	1.16689E-02	3.46944E+01	3.22044E+05	2.96149E+00	1.15147E+02	1.07138E+01	1.09162E-01	4.41139E+06	2.66234E+01	4.65234E-01	6.42001E+06	7.84685E+02	1.22753E+05	5.75853E+00	1.90942E+02	9.74171E+00	4.88446E+06];



B = model.p_nominal;
	
for i=1:length(model.p_estimate)
	B(:,(model.p_estimate(i)))= p(:,i);
end

samples=1000;

for pp=1:samples
	disp(pp);        
	xInit = samplingScheme(model.x_init,'uniform',0.05);
	tf = Gen_Graph(B,xInit,nTime,dt);
	AB{pp} = transpose(reshape(tf,length(model.x_init),nTime));
end

t= linspace(0,nTime*dt,nTime);
save('egfngf_paramest');

%%
load('egfngf_paramest');
figure; hold on;
idx = [4 6 8 27 25 10 21];
chooseidx = [2 3 5 7];
xidx = idx(chooseidx);
nplot = length(xidx);
for i=1:nplot
	subplot(1,nplot,i);		
	for pp=1:samples
		axis([0 nTime 0 1.2*max(max(AB{pp}(:,xidx(i))),max(numData(:,chooseidx(i)+1)))]);
		plot(t, AB{pp}(:,xidx(i)),'Color',[0.6,0.6,0.6]);
		hold on;
	end
	
	title(model.x_names{xidx(i)},'FontSize',14);

	errorbar(numData(:,1),numData(:,chooseidx(i)+1),numData(:,chooseidx(i)+1)*0.1,'ks','Color',[0 0 0],'MarkerEdgeColor','k',...
	 'MarkerFaceColor','w','MarkerSize',5);
end
%saveas(h,strcat('../../models/egf-ngf/test/sres-1-train',model.x_names{str2num(textData{i}(2:end))}),'pdf');
%close
%%
idx2 = [19 23];
chooseidx2 = [2];
xidx2 = idx2(chooseidx2);
nplot2 = length(xidx2);
figure; hold on;

for i=1:nplot2
 subplot(1,nplot2,i);
 for pp=1:samples
	axis([0 nTime 0 1.2*max(max(AB{pp}(:,xidx2(i))),max(numData2(:,chooseidx2(i)+1)))]);
	plot(t, AB{pp}(:,xidx2(i)),'Color',[0.6,0.6,0.6]);
	hold on;
 end

 %xlabel('Model time','FontSize',14);
 title(strcat(model.x_names{xidx2(i)}),'FontSize',14);

 errorbar(numData2(:,1),numData2(:,chooseidx2(i)+1),numData2(:,chooseidx2(i)+1)*0.1,'ks','MarkerEdgeColor','k',...
	 'MarkerFaceColor','w','MarkerSize',5);
	 hold on;

end