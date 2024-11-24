#! /bin/bash

usr="$(id -u -n)"
initial_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
run_nums="53285"
# run_nums="41989 41992"
# run_nums="53341 53376 53370 53387"
# run_nums="52911"
# file_format="/sphenix/lustre01/sphnxpro/commissioning/slurp/tpcbeam/run_%08d_%08d/*%08d*"
file_format="/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/new_2024p002/run_%08d_%08d/*%08d*"

cd $initial_dir
list="segments.list"

rm -rf "$list"
for run_num in $run_nums; do
# for run_num in "$@"; do
	file_pattern=$(printf $file_format $(( (run_num / 100) * 100 )) $(( (run_num / 100 + 1) * 100 )) $run_num)
	echo "$file_pattern"
	ls -1 $file_pattern >> $list 2> /dev/null
done

line_count=$(wc -l $list | awk '{print $1}')

job="job.job"
cat << EOF > $job
universe           = vanilla
executable         = wrapper.sh
arguments          = \$(Process) $initial_dir/$list
initialdir         = $initial_dir

notification       = Never

output             = $initial_dir/out/out_\$(Process).out
error              = $initial_dir/out/out_\$(Process).out
log                = /tmp/${usr}_alignment_\$(Process).log

initialdir         = $initial_dir
request_memory     = 8192MB
PeriodicHold       = (NumJobStarts >= 1 && JobStatus == 1)
concurrency_limits = CONCURRENCY_LIMIT_DEFAULT:100

queue $line_count
EOF

echo "submitting job file $job"
condor_submit $job

