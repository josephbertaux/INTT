#! /bin/bash

export USER="$(id -u -n)"
export LOGNAME="${USER}"
export HOME="/sphenix/u/${LOGNAME}"

initial_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
# export MYINSTALL=""
export MYINSTALL="/sphenix/user/jbertaux/MYINSTALL"

if [ $# -ne 2 ]; then
cat << EOF

	usage:
		$0 [line number] [file list]

EOF
	exit 0
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
if [ -n "$MYINSTALL" ] && [ -d "$MYINSTALL" ]; then
	source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
fi

cd $initial_dir

root -q -b "Fun4All_FieldOnAllTrackers.C($1, \"$2\")"
# gdb -ex run --args root.exe -q -b "Fun4All_FieldOnAllTrackers.C($1, \"segments.list\")"

rv=$?
exit $rv

