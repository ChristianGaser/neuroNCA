# Personal Makefile variables
#
# $Id$

OLDVERSION="0.9"
NEWVERSION="1.0"
REVISION=`git rev-list --count HEAD`
DATE=`git log --date short |grep "Date:"|head -1|cut -f2 -d':'|sed -e s'/ //g'`

ZIPFOLDER=/Users/gaser/matlab/NCA

TARGET=/Users/gaser/spm/spm12/toolbox/NCA
TARGET2=/Volumes/UltraMax/spm12/toolbox/NCA
TARGET3=paris.biomag.uni-jena.de:/home/gaser/spm12/toolbox/NCA

STARGET_HOST=141.35.69.218
STARGET_HTDOCS=${STARGET_HOST}:/volume1/web/
STARGET_FOLDER=/volume1/web/NCA

STARGET=${STARGET_HOST}:${STARGET_FOLDER}

MATLAB_FILES=*NCA*.m snpm_P_FDR.m

FILES=${MATLAB_FILES} ${C_FILES} ${MISC_FILES}

ZIPFILE=neurNCA_r${REVISION}.zip

install:
	-@echo install
	-@test ! -d ${TARGET} || rm -rf ${TARGET}
	-@mkdir ${TARGET}
	-@cp -R ${FILES} ${TARGET}

install2:
	-@echo install2
	-@test ! -d ${TARGET2} || rm -rf ${TARGET2}
	-@mkdir ${TARGET2}
	-@cp -R ${FILES} ${TARGET2}

install3:
	-@echo install3
	-@scp -r ${FILES} ${TARGET3}/

help:
	-@echo Available commands:
	-@echo install install2 install3 zip scp doc update

update:
	-@git fetch
	-@echo '% NCA Toolbox' > Contents.m
	-@echo '% Version ' ${REVISION} ' (version '${NEWVERSION}')' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m
	-@perl -p -i -e "s/${OLDVERSION}/${NEWVERSION}/g" spm_NCA.m
	-@echo '% NCA Toolbox' > INSTALL.txt
	-@echo '% Version ' ${REVISION} ${NEWVERSION} ${DATE} >> INSTALL.txt
	-@cat INSTALL_info.txt >> INSTALL.txt

zip: update
	-@echo zip
	-@test ! -d NCA || rm -r NCA
	-@mkdir NCA
	-@cp -rp ${FILES} NCA
	-@bash update_revision.sh
	-@zip ${ZIPFOLDER}/${ZIPFILE} -rm NCA

scp: zip
	-@echo scp to http://${STARGET_HOST}/NCA/${ZIPFILE}
	-@scp -O -P 2222 CHANGES.txt ${ZIPFOLDER}/${ZIPFILE} ${STARGET}
	-@bash -c "ssh -p 2222 ${STARGET_HOST} ln -fs ${STARGET_FOLDER}/${ZIPFILE} ${STARGET_FOLDER}/NCA_latest.zip"
	
