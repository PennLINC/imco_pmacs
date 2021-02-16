#!/bin/csh
######################
##### 02/16/2021 #####
######################

#pre: needs a csv with final sample, first column bblid, second column date x scanid
#post: average lh and rh maps for alff and cbf individually
#uses: When we moved alff and cbf surface maps over from chead, it included all surface projections, including ones not in our analysis. To make it easy to make a group average, this script
 #1) makes a directory within the alff (or cbf) directories,
 #2) symbolically links to the files in the n831 group only, and
 #3) runs the freesurfer command mri_concat with --mean as an option, to create a mean lh and rh map.
 #4) These are then converted to asc.
 #5) the asc files are saved in the results section, with asc swapped for csv suffix for use in the matlab visualization script
 #6) coupling maps are also converted to ascii and moved to mean maps, also with .csv suffixes
#dependenccies: freesurfer (pmacs uses version 7.1.1
#makes the average files for display in imco paper

#make results destination
mkdir /project/imco/baller/results/mean_maps

#first, do with alff
cd /project/imco/surfaceMaps/alff_from_chead/
mkdir n831_maps
cd n831_maps
rm *
set subj_num = 0
foreach subj (`cat /project/imco/baller/subjectLists/n831_alff_cbf_finalSample.csv`)
    if ($subj_num > 0) then
	set bblid = `echo $subj | cut -f1 -d ','`
	set scanid = `echo $subj | cut -f2 -d ','`
	ln -s ../${bblid}_${scanid}*.mgh . 
    endif
    @ subj_num+=1
end

# make average maps
mri_concat --i *lh*.mgh --o lh_alff_mean.mgh --mean
mri_concat --i *rh*.mgh --o rh_alff_mean.mgh --mean

#convert to asc
mri_convert --ascii lh_alff_mean.mgh lh_alff_mean.asc
mri_convert --ascii rh_alff_mean.mgh rh_alff_mean.asc

#copy to results directory
cp lh_alff_mean.asc /project/imco/baller/results/mean_maps/lh_alff_mean.csv
cp rh_alff_mean.asc /project/imco/baller/results/mean_maps/rh_alff_mean.csv

#### now do same with cbf###
cd /project/imco/surfaceMaps/cbf_from_chead/
mkdir n831_maps
cd n831_maps
rm *
set subj_num = 0
foreach subj (`cat /project/imco/baller/subjectLists/n831_alff_cbf_finalSample.csv`)
    if ($subj_num > 0) then
	set bblid = `echo $subj | cut -f1 -d ','`
	set scanid = `echo $subj | cut -f2 -d ','`
	echo $bblid $scanid
	ln -s ../${bblid}_${scanid}*.mgh . 
    endif
    @ subj_num+=1
end
mri_concat *lh*.mgh --o lh_cbf_mean.mgh --mean
mri_concat *rh*.mgh --o rh_cbf_mean.mgh --mean

#convert to asc
mri_convert --ascii lh_cbf_mean.mgh lh_cbf_mean.asc
mri_convert --ascii rh_cbf_mean.mgh rh_cbf_mean.asc

#copy to results directory
cp lh_cbf_mean.asc /project/imco/baller/results/mean_maps/lh_cbf_mean.csv
cp rh_cbf_mean.asc /project/imco/baller/results/mean_maps/rh_cbf_mean.csv

#copy coupling average maps into mean maps for good measure
cd /project/imco/couplingSurfaceMaps
mri_convert --ascii n831_lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.mgh n831_lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc
mri_convert --ascii n831_rh.coupling_coef_alff_cbf.fwhm15.fsaverage5.mgh n831_rh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc

cp n831_lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc /project/imco/baller/results/mean_maps/n831_lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.csv
cp n831_rh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc /project/imco/baller/results/mean_maps/n831_rh.coupling_coef_alff_cbf.fwhm15.fsaverage5.csv

#reset directory to be back at the beginning
cd /project/imco/baller/scripts
