#! /bin/bash

export USER="$(id -u -n)"
export LOGNAME="${USER}"
export HOME="/sphenix/u/${USER}"
export MYINSTALL="/sphenix/user/${USER}/MYINSTALL"

HELP(){
cat << EOF
	
	usage:
		$0 [/full/path/to/macro.C] [subprocess id] [number of events]
	Wrapper shell script to run either Fun4All macro

EOF
}

if [ ! -f "$1" ] || ! [[ $2 =~ ^[0-9]+$ ]] || ! [[ $3 =~ ^[0-9]+$ ]]; then
	HELP
	exit 0
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
if [ -n "${MYINSTALL}" ] && [ -d "${MYINSTALL}" ]; then
	source /opt/sphenix/core/bin/setup_local.sh ${MYINSTALL}
fi

USE_CONDOR=0
if [ -n "${_CONDOR_SCRATCH_DIR}" ] && [ -d "${_CONDOR_SCRATCH_DIR}" ]; then
	cd ${_CONDOR_SCRATCH_DIR}
	rsync -av $(dirname $1) .
	USE_CONDOR=1
else
	cd $(dirname $1)
fi

cat << EOF
	
	Running in directory
		$(pwd)
	ls:
		$(ls -l)

EOF

echo root -q -b "$(basename $1)(\"$2\", $3, ${USE_CONDOR})"
# root -q -b "$(basename $1(\"$2\", $3, ${USE_CONDOR}))"
# gdb -ex run --args root.exe -q -b "$(basename $1)(\"$2\", $3, ${USE_CONDOR})"
RV=$?


if [ $USE_CONDOR -eq 1 ]; then
	echo "asdf"
	[ -d "dataloader" ] && cp -r "$(dirname $1)/."
	[ -d "factory" ] && cp -r "$(dirname $1)/."
	cp "*.root" "/sphenix/tg/tg01/hf/${USER}/."
fi

exit $RV

