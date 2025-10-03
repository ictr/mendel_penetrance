# 2020.01.17.1107: Per Chris, installing mendel
# but have to build gcc with fortran etc. first

setup1st:
	@echo "Make sure you do the following first:"
	@echo "module load gcc/9.2.0"
	@echo "Next: make go"

STEP1_BINARY = step1/step1_single


step1-compile: $(STEP1_BINARY)

$(STEP1_BINARY): step1/lfs_single.f
	cd step1 && gfortran -O3 mendela-batch.f lfs_single.f -o step1_single

step1-example: $(STEP1_BINARY)
	cd data && ../$(STEP1_BINARY) && mv single_out.dat data_step1.dat && cat data_step1.dat | grep -v IERROR | grep -vx ''

step1-crc: $(STEP1_BINARY)
	cd colorectal_raw && ../$(STEP1_BINARY) && mv single_out.dat crc_step1.dat && cat crc_step1.dat | grep -v IERROR | grep -vx ''

step1-breast: $(STEP1_BINARY)
	cd breast_raw && ../$(STEP1_BINARY) && mv single_out.dat breast_step1.dat && cat breast_step1.dat | grep -v IERROR | grep -vx ''

STEP2_CRC_BINARY = step2/step2_crc
STEP2_BREAST_BINARY = step2/step2_breast

$(STEP2_CRC_BINARY): step2/lfs_crc.f
	cd step2 && gfortran -O3 mendela-batch.f lfs_crc.f -o step2_crc

$(STEP2_BREAST_BINARY): step2/lfs_breast.f
	cd step2 && gfortran -O3 mendela-batch.f lfs_breast.f -o step2_breast

step2-crc: $(STEP2_CRC_BINARY)
	cd colorectal_raw && ../$(STEP2_CRC_BINARY) && mv single_out.dat crc_step2.dat && cat crc_step2.dat | grep -v IERROR | grep -vx ''

step2-breast: $(STEP2_BREAST_BINARY)
	cd breast_raw && ../$(STEP2_BREAST_BINARY) && mv single_out.dat breast_step2.dat && cat breast_step2.dat | grep -v IERROR | grep -vx ''


STEP3_CRC_BINARY = step3/step3_crc
STEP3_BREAST_BINARY = step3/step3_breast
STEP3_CRC_SRC = step3/lfs_crc.f
STEP3_BREAST_SRC = step3/lfs_breast.f
STEP3_TEMPLATE = step3/single_template.f
STEP3_SEX_SPECIFIC_TEMPLATE = step3/single_sex_specific_template.f
STEP3_CRC_PARAMS = step3/step3_crc_params.json
STEP3_BREAST_PARAMS = step3/step3_breast_params.json

STEP4_TEMPLATE = step3/single_template.f
STEP4_BREAST_PARAMS = step4/step4_breast_params.json
STEP4_BREAST_BINARY = step4/step4_breast

${STEP3_CRC_SRC}: ${STEP3_TEMPLATE} ${STEP3_CRC_PARAMS}
	@echo "Generating step3/lfs_crc.f from template"
	cat ${STEP3_CRC_PARAMS} | jinja2 ${STEP3_TEMPLATE} --format=json  > ${STEP3_CRC_SRC}

${STEP3_BREAST_SRC}: ${STEP3_TEMPLATE} ${STEP3_BREAST_PARAMS}
	@echo "Generating step3/lfs_breast.f from template"
	cat ${STEP3_BREAST_PARAMS} | jinja2 ${STEP3_TEMPLATE} --format=json  > ${STEP3_BREAST_SRC}


$(STEP3_CRC_BINARY): ${STEP3_CRC_SRC}
	cd step3 && gfortran -O3 mendela-batch.f lfs_crc.f -o step3_crc

$(STEP3_BREAST_BINARY): ${STEP3_BREAST_SRC}
	cd step3 && gfortran -O3 mendela-batch.f lfs_breast.f -o step3_breast


step3-crc: $(STEP3_CRC_BINARY)
	cd colorectal_raw && ../$(STEP3_CRC_BINARY) && mv single_out.dat crc_step3.dat && cat crc_step3.dat | grep -v IERROR | grep -vx ''

step3-breast: $(STEP3_BREAST_BINARY)
	cd breast_raw && ../$(STEP3_BREAST_BINARY) > step3_breast_debug.txt
	cd breast_raw && mv single_out.dat breast_step3.dat && cat breast_step3.dat | grep -v IERROR | grep -vx ''

${STEP4_BREAST_SRC}: ${STEP4_TEMPLATE} ${STEP4_BREAST_PARAMS}
	@echo "Generating step4/lfs_breast.f from template"
	cat ${STEP4_BREAST_PARAMS} | jinja2 ${STEP4_TEMPLATE} --format=json  > ${STEP4_BREAST_SRC}

$(STEP4_BREAST_BINARY): ${STEP4_BREAST_SRC}
	cd step4 && gfortran -O3 mendela-batch.f lfs_breast.f -o step4_breast

step4-breast: $(STEP4_BREAST_BINARY)
	cd breast_raw && ../$(STEP4_BREAST_BINARY) > step4_breast_debug.txt
	cd breast_raw && mv single_out.dat breast_step3.dat && cat breast_step3.dat | grep -v IERROR | grep -vx ''


clean:
	\rm -rf single poly
