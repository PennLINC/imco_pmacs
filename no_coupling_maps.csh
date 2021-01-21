#!/bin/csh
set cnt = 1
echo -n "" >! /project/imco/baller/subjectLists/n831_no_freesurfer_dir
foreach x (`cat /project/imco/baller/subjectLists/n831_imageOrder.csv`)
	if ($cnt > 1) then
	   set bblid = `echo $x | cut -f1 -d ','`
 	   set datexscanid = `echo $x | cut -f2 -d ','`
	   echo $bblid/$datexscanid
	   if (`ls /project/imco/surfaces/$bblid/$datexscanid/surf/lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc | wc -l` > 0) then
	    echo "$bblid/$datexscanid exists"
	   else
	    echo $bblid/$datexscanid" does not exist"
	    echo $bblid/$datexscanid >> /project/imco/baller/subjectLists/n831_no_freesurfer_dir
	   endif
	endif
	@ cnt += 1 
end
