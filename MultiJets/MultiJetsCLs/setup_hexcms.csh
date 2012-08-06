setenv VO_CMS_SW_DIR /cms/base/cmssoft
setenv SCRAM_ARCH slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.csh

source /osg/apps/osg/setup.csh

setenv PATH /home/cdfcaf/condor/dist/bin:${PATH}

setenv ROOTSYS /cms/data24/dhidas/CernRoot/root_2011.09.17_HEAD/root
setenv LD_LIBRARY_PATH $ROOTSYS/lib
set path = ( $path $ROOTSYS/bin )


#cd /users/h2/dhidas/CMSSW_4_2_5/src
#eval `scramv1 runtime -csh`
#cd -


