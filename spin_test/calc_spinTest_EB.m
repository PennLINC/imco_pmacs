% addpaths
addpath(genpath('/project/imco/baller/scripts/spin_test/'));
ProjectFolder = '/project/imco/baller/results/coupling_accuracy/';
%Variability_Visualize_Folder = [WorkingFolder '/Variability_Visualize'];
% set outdir
outdir= '/project/imco/baller/results/coupling_accuracy/spin_test_results/';
% load in trinarized masks, 1 for FDR and in yeo masks, 0 for not fdr, -1 for in medial wall. 
%surfML = fileread([ProjectFolder, '/lh_lm_sex_t_fdr05_Yeo7_1_0_-1_spin.csv']);
%surfMR = fileread([ProjectFolder, '/rh_lm_sex_t_fdr05_Yeo7_1_0_-1_spin.csv']);

%surfML = '/cbica/projects/pinesParcels/data/H_SNR_masks/lh.Mask_SNR.label';
%mwIndVec_l = read_medial_wall_label(surfML);
%surfMR = '/cbica/projects/pinesParcels/data/H_SNR_masks/rh.Mask_SNR.label';
%mwIndVec_r = read_medial_wall_label(surfMR);
% for each scale
%for K=2:30
	% report K
%	disp(K)
	% set output file name
outFn=strcat([outdir, 'spin_results_yeo_1_0_-1_output']);
	% load in MAD values
 %       MADFileP = ['/gpfs/fs001/cbica/projects/pinesParcels/data/SingleParcellation/SingleAtlas_Analysis/Variability_Visualize/VariabilityLoading_Median_' num2str(K) 'SystemMean.mat'];
  %      MADFile=load(MADFileP);
%	mad_lh=MADFile.VariabilityLoading_Median_KSystemMean_lh
%	mad_rh=MADFile.VariabilityLoading_Median_KSystemMean_rh
	
	% set mask ROI values to -1 for spin test to catch em as invalid
%	mad_lh(mwIndVec_l)=100;
%	mad_rh(mwIndVec_r)=100;
	% write them out as a transposed csv for spin test to deal with		
%	writetable(table(mad_lh'),[outdir num2str(K) '_L.csv'],'WriteVariableNames',0);
%	writetable(table(mad_rh'),[outdir num2str(K),'_R.csv'],'WriteVariableNames',0);
	% create permutations, save out to outFn
SpinPermuFS([ProjectFolder, '/lh_lm_sex_t_fdr05_Yeo7_1_0_-1_spin.csv'],[ProjectFolder, '/rh_lm_sex_t_fdr05_Yeo7_1_0_-1_spin.csv'], 1000, outFn);
%end
