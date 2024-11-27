#!/bin/bash

USR="$(id -u -n)"
NEVT=2000

# The location of this shell script
PWD=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

mkdir -p ${PWD}/job
mkdir -p ${PWD}/out

# used in naming the job file
# FILE="${PWD}/job/$(basename ${EXE} .sh)_${ARGS}.job"
FILE="${PWD}/job/sig_jobs.job"
cat << EOF > ${FILE}
universe           = vanilla
executable         = ${PWD}/scripts/gen.sh
arguments          = ${PWD}/macro/Fun4All_HF.C \$(process) ${NEVT}

notification       = Never

output             = ${PWD}/out/gen_sig_\$(process)_${NEVT}.out
error              = ${PWD}/out/gen_sig_\$(process)_${NEVT}.out
log                = /tmp/${USR}_gen_sig_\$(process)_${NEVT}.log

request_memory     = 8192MB
PeriodicHold       = (NumJobStarts >= 1 && JobStatus == 1)
concurrency_limits = CONCURRENCY_LIMIT_DEFAULT:100

queue 3
EOF

echo "submitting job file ${FILE}"
condor_submit ${FILE}

exit 0
