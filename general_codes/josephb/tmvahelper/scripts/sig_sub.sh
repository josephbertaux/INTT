#!/bin/bash

NUM_EVT=2000
NUM_JOB=1

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source tmva_setup.sh

mkdir -p ${TMVA_SOURCE_DIR}/job
mkdir -p ${TMVA_SOURCE_DIR}/out

# used in naming the job file
FILE="${TMVA_SOURCE_DIR}/job/sig_jobs.job"
cat << EOF > ${FILE}
universe           = vanilla
executable         = ${TMVA_SOURCE_DIR}/scripts/gen.sh
arguments          = ${TMVA_SOURCE_DIR}/macro/Fun4All_HF.C \$(process) ${NUM_EVT}

notification       = Never

output             = ${TMVA_SOURCE_DIR}/out/gen_sig_\$(process)_${NUM_EVT}.out
error              = ${TMVA_SOURCE_DIR}/out/gen_sig_\$(process)_${NUM_EVT}.out
log                = /tmp/${USER}_gen_sig_\$(process)_${NUM_EVT}.log

request_memory     = 8192MB
PeriodicHold       = (NumJobStarts >= 1 && JobStatus == 1)
concurrency_limits = CONCURRENCY_LIMIT_DEFAULT:100

queue ${NUM_JOB}
EOF

echo "submitting job file ${FILE}"
condor_submit ${FILE}

exit 0
