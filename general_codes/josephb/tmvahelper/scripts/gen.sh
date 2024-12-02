#! /bin/bash

SHOW(){
cat << EOF
Running in directory:
	$(pwd)

Directory contents:
$(ls -la)

EOF
}

if [ ! -f "$1" ] || ! [[ $2 =~ ^[0-9]+$ ]] || ! [[ $3 =~ ^[0-9]+$ ]]; then
cat << EOF
	
	usage:
		$0 [macro.C] [subprocess id] [number of events]
	Wrapper shell script to run either Fun4All macro which generates training data
	(macro/Fun4All_HF.C for signal, macro/Fun4All_MB.C for background)

	Must be run as a condor job

EOF
	exit 0
fi

if [ -z "${_CONDOR_SCRATCH_DIR}" ] || ! [ -d "${_CONDOR_SCRATCH_DIR}" ]; then
cat << EOF

	Job must run under condor

EOF
	exit 0
fi

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source tmva_setup.sh
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
if [ -n "${MYINSTALL}" ] && [ -d "${MYINSTALL}" ]; then
	source /opt/sphenix/core/bin/setup_local.sh ${MYINSTALL}
fi

cd ${_CONDOR_SCRATCH_DIR}
rsync -av ${TMVA_SOURCE_DIR}/macro .
SHOW

root -q -b "$(basename $1)(\"$2\", $3)"
# gdb -ex run --args root.exe -q -b "$(basename $1)(\"$2\", $3)"
RV=$?

SHOW
cp *.root ${TMVA_DATA_DIR}/.

echo "$0" done
exit $RV

