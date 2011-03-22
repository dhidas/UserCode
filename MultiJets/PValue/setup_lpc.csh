source /uscmst1/prod/sw/cms/cshrc prod
source /uscmst1/prod/grid/gLite_SL5.csh
cmscvsroot CMSSW

cd /uscms/home/dhidas/CMSSW_3_10_0/src
eval `scramv1 runtime -csh`
cd -

set arch = slc5_ia32_gcc434
set roofit_release = 5.27.06

setenv CMS_PATH /uscmst1/prod/sw/cms

setenv PYTHONDIR /uscmst1/prod/sw/cms/slc5_ia32_gcc434/external/python/2.6.4-cms8

setenv PATH ${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.27.06/bin:$PATH
setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.27.06
setenv PYTHONPATH ${ROOTSYS}/lib

setenv LD_LIBRARY_PATH ${PYTHONDIR}/lib:${CMS_PATH}/$arch/external/gcc/4.3.4/lib:${ROOTSYS}:${ROOTSYS}/lib:${LD_LIBRARY_PATH}

