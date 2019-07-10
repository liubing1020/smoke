clear variables
addpath ode/
addpath util/
addpath analysis/
addpath ../models/repressilator/

cd ode
compile('repress','graph')
cd ..
model = repressilator_pathway();

%% Define discretization intervals
% State discretization
% Time between model checking steps (timse discretization)
dt = 0.5;
% Number of discrete time steps
nTime = 30;

%sres-normal
%p=[85.2018	86.6885	88.5489	190.4744	159.4181	189.6430	12.4335	10.4924	13.4802];
%ga-normal
%p=[73.0524	70.8285	60.5383	158.8891	149.7040	127.9257	9.3535	7.6987	10.7029];
%ga-pval
%p=[40.5817	12.3458	58.4368	106.6602	21.6903	117.5679	11.1524	8.1342	8.0321];
%sres-pval
%p=[85.1120	72.4767	96.4276	182.2802	158.3923	199.8465	10.7862	8.6189	8.4370];
%0.8-p
%p=[7.14383E+01	6.95836E+01	7.26164E+01	1.52564E+02	1.39507E+02	1.54691E+02	1.19479E+01	8.76358E+00	1.04238E+01];
%0.8
%p=[8.12189E+01	5.19553E+01	7.55776E+01	1.89710E+02	8.80473E+01	1.63956E+02	1.08700E+01	8.12559E+00	1.19910E+01];
%0.9
%p=[8.61548E+01	6.89089E+01	7.31270E+01	1.78593E+02	1.30040E+02	1.56908E+02	1.17314E+01	1.21534E+01	1.44455E+01];
%0.9-p
p=[8.00087E+01	9.20495E+01	5.61409E+01	1.68510E+02	1.76316E+02	1.40932E+02	6.88332E+00	7.52111E+00	1.26742E+01];
%GA
%0.8-p
%p=[2.34568E+01	7.78209E+01	8.19246E+01	4.76258E+01	1.45943E+02	1.62647E+02	1.23930E+01	1.23638E+01	1.22635E+01];
%0.8
%p=[3.36533E+01	5.92057E+01	6.25271E+01	7.87844E+01	1.00586E+02	1.28899E+02	1.16558E+01	6.61712E+00	1.08631E+01];
%0.9-p
%p=[6.13435E+01	7.62342E+01	4.83178E+01	1.31852E+02	1.44284E+02	1.03991E+02	1.02796E+01	9.42117E+00	1.13869E+01];
%0.9
%p=[2.05680E+01	8.91487E+01	8.11095E+01	4.90906E+01	1.58555E+02	1.67449E+02	1.28105E+01	7.02243E+00	1.50332E+01];

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
a=[1 2 3 4 5 6];

%%

h=figure();
p=1;
xnames = {'m_1','m_2','m_3','p_1','p_2','p_3'};
for k=1:length(model.x_init)
	if find(k==a)
		subplot(2,3,p);	
		p=p+1;
		for pp=1:samples
			plot(t,AB{pp}(:,k),'LineWidth',2);
			hold on;
		end
		xlabel('Model time','FontSize',14);
		ylabel(xnames{k},'FontSize',14);
		set(gca,'FontSize',12)
	end
end