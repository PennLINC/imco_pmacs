#!/bin/csh
cd /project/imco/surfaceMaps/alff_from_chead/
foreach x (`ls *.mgh*`)
	set asc_name = `echo $x | perl -pe 's/mgh/asc/'`
	echo "$x : $asc_name"
	mri_convert --ascii+crsf $x $asc_name
end

cd /project/imco/surfaceMaps/cbf_from_chead
foreach x (`ls *.mgh*`)
    set asc_name = `echo $x | perl -pe 's/mgh/asc/'`
    echo "$x : $asc_name"
    mri_convert --ascii+crsf $x $asc_name
end
cd /project/imco/baller/scripts/
