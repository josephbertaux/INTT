#!/bin/bash

USR="$(id -u -n)"
PWD=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# The location of this shell script

EXE="${PWD}/shell.sh"
DIR="${PWD}"
MEM="4096MB"

show_help() {
cat << EOF
	usage:
		./$0 [run number] [num events (optional)] [run type (optional)]
	submits a condor job to do a calibration for [run number]
	the output will be in the local directory
EOF
}

if [[ $# -lt 1 ]]; then
	show_help
	exit 1
fi

if [[ $1 == "-h" ]]; then
	show_help
	exit 0
fi

RUN_NUM="$1"

[[ -d ${DIR}/job ]] || mkdir -p ${DIR}/job
[[ -d ${DIR}/out ]] || mkdir -p ${DIR}/out
[[ -d ${DIR}/err ]] || mkdir -p ${DIR}/err
[[ -d ${DIR}/log ]] || mkdir -p ${DIR}/log

# rm -f ${DIR}/out/*
# rm -f ${DIR}/err/*

# underscore-deliminated list of command line args
IFS="_"
ARGS="$*"
IFS=" "

# used in naming the job file
# FILE="${DIR}/job/$(basename ${EXE} .sh)_${ARGS}.job"
FILE="${DIR}/job/job_${ARGS}.job"
cat << EOF > ${FILE}
universe        = vanilla
executable      = ${EXE}
arguments       = \$(Process) $*

notification    = Never

output          = ${DIR}/out/out_\$(Process)_${ARGS}.out
error           = ${DIR}/err/err_\$(Process)_${ARGS}.err
log             = ${DIR}/log/log_\$(Process)_${ARGS}.log

initialdir      = ${PWD}
request_memory  = ${MEM}
periodichold    = (NumJobStarts >= 1 && JobStatus == 1)

queue 8
EOF

echo "submitting job file ${FILE}"
condor_submit ${FILE}

exit 0
