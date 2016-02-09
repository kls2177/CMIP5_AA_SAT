#!/bin/csh
set echo

#EXTRACT_2DVAR_CMIP5: This is a wrapper script for extracting 2D CMIP5 fields from local CMIP5 repository.
#Script extracts the field and then passes it to python to interpolate onto 2x2 lat-lon grid
    
#User defined fields
set savdir = /home/ksmith/climdata/CMIP5/mon/piControl
set pwdsed = `echo $savdir | sed 's_/_\\/_g'`
set exedir = /home/ksmith/Insight
set varname = ("tas")
#===============================================
cd $savdir

###loop over models
set nmod = `ls | wc -l`
echo $nmod
foreach imod (`seq 1 $nmod`)
@ modline = $imod + 1
set modname = `ls -l | sed -n $modline'p' | awk '{print $9}'`
cd $modname

###loop over ensemble members
set nmem = `ls | wc -l`
echo $nmem
foreach imem (`seq 1 $nmem`)
@ memline = $imem + 1
set memname = `ls -l | sed -n $memline'p' | awk '{print $9}'`
cd $memname

set ct = `\ls -l $savdir/$modname/$memname/$varname.*.nc | wc -l`
if ( $ct == 0 ) then  ## If there are no files to process, then go back to root
   cd $savdir/$modname
else
# call python module
sed s/XMODELX/${modname}/ $exedir/AA_SAT_CMIP5_pictrl_diags.py > tmp.py
sed s/XMEMX/${memname}/ tmp.py > tmp2.py
sed s/XVARX/${varname}/ tmp2.py > tmp3.py
sed s/XDIRX/${pwdsed}/ tmp3.py > tmp4.py
nohup /home/ksmith/anaconda3/bin/python3 < tmp4.py > run.out &
wait	
rm -rf tmp*.py
cd $savdir/$modname
endif # if statment
end # loop over ensemble members
cd $savdir
end # loop over models
