setenv VO_CMS_SW_DIR /cms/base/cmssoft
setenv SCRAM_ARCH slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.csh

source /osg/apps/osg/setup.csh

setenv PATH /home/cdfcaf/condor/dist/bin:${PATH}

cd /users/h2/dhidas/CMSSW_4_2_5/src
eval `scramv1 runtime -csh`
setenv PATH /cms/data25/krose/RootInstall/root/bin:${PATH}
setenv LD_LIBRARY_PATH /cms/data25/krose/RootInstall/root/lib:${LD_LIBRARY_PATH}
cd -


