% addpaths
%addpath(genpath('/project/imco/baller/scripts/spin_test/functions/'));
% set outdir
outdir='/project/imco/baller/results/coupling_accuracy/spin_test_results/';
% set paths
ProjectFolder = '/project/imco/baller/results/coupling_accuracy/spin_test_results';
% get trinarized vectors, 1 for in yeo network, 0 for not, -1 for in medial wall
%vectl = fileread([ProjectFolder, '/lh_lm_sex_t_fdr05_Yeo7_1_0_-1_spin.csv']);
%vectr = fileread([ProjectFolder, '/rh_lm_sex_t_fdr05_Yeo7_1_0_-1_spin.csv']);

% initialize permutation house for proportions for 1000 spins across scales, +1 row for real correlation
permHouse=zeros(1001,7);
Index_l=[1:10242];
Index_r=[1:10242];
% for each scale, get disitribution of spatial correlations with PG1
for K=2:8 % parcels 1-7, have no idea if this is right. 
%	disp(K)
	% get MAD to test
	%MADFile = strcat([outdir, 'spin_results_yeo_1_0_-1_output.mat']);
%['/gpfs/fs001/cbica/projects/pinesParcels/data/SingleParcellation/SingleAtlas_Analysis/Variability_Visualize/VariabilityLoading_Median_' num2str(K) 'SystemMean.mat'];
	%initmat=load(MADFile);
	%MAD_atK=initmat.VariabilityLoading_Median_KSystemMean_NoMedialWall;
	%% get real correlation
	%realrho=corr(MAD_atK',pg1,'type','spearman','rows','complete');
	%permHouse(1,(K-1))=realrho;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% get permuted correlations
	permutFile=strcat([outdir, 'spin_results_yeo_1_0_-1_output.mat']);
	permuts=load(permutFile);
	permutsL=permuts.bigrotl(:,Index_l);
	permutsR=permuts.bigrotr(:,Index_r);	
	% change 100 (markers of invalid vertices) to NA
	permutsL(permutsL==-1)=NaN;
	permutsR(permutsR==-1)=NaN;

	%load in all indexes DA, FP, etc

	
	% for each permutation
	for P=1:1000
		permutVals=[permutsL(P,:) permutsR(P,:)];

		% figure out network assignment for each
		total_vals = (DMN)
		medial_wall = (#-1s in that parcel)
		% remove the medial wall stuff
		denominator = (total_vals - medial_wall_vals)
		% take the proportion
		permutProp = permutVals/denominator
		% store
		permHouse((1+P),(K-1))=permProp;
%		permrho=corr(permutVals',MAD_atK','type','spearman', 'rows','complete');
		
	end
end
% write out distribution, R friendly format
writetable(array2table(permHouse),strcat(outdir,'SpinTestDistrs_MAD_PG1.csv'),'Delimiter',',','QuoteStrings',true);
