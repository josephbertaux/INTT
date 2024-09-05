#!/bin/bash

# MYINSTALL=""
MYINSTALL="/sphenix/user/jbertaux/MYINSTALL"

# One-liner to get pwd
PWD=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

INTT_FORMAT="/sphenix/lustre01/sphnxpro/physics/INTT/%s/%s_intt%d-%08d-*"
LIST_FORMAT="${PWD}/lst/run_%08d_intt%01d.list"
DATA_FORMAT="${PWD}/biz/run_%08d_intt%01d.root"

RUN_NUM=""         # required
NUM_EVT="10000000" # default if argument 2 is empty
RUN_TYPE="physics" # default if argument 3 is empty

show_help()
{
cat << EOF

	usage:
		$0 [which intt] [run number] [num events (optional)] [run type (optional)]
	The felix server this should run for is specified by [which intt]
	The default for [num events] is ${NUM_EVT}
	The default for [run type] is ${RUN_TYPE}
	For older runs, you may need to set [run type] as "beam" (and specify [num events])

EOF
}

# Arguments
if [[ $# -lt 1 || $1 == "-h" || $1 == "--help" ]]; then
	show_help
	exit 0
fi

WHICH_INTT="$1"

RUN_NUM="$2"

if [[ -n "$3" ]]; then
	NUM_EVT="$3"
fi

if [[ -n "$4" ]]; then
	RUN_TYPE="$4"
fi

# Custom MYINSTALL for developement/debugging purposes
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
if [ -n "${MYINSTALL}" ] && [ -d "${MYINSTALL}" ]; then
	source /opt/sphenix/core/bin/setup_local.sh ${MYINSTALL}
fi

printf -v LIST ${LIST_FORMAT} ${RUN_NUM} ${WHICH_INTT}
printf -v FILE ${INTT_FORMAT} ${RUN_TYPE} ${RUN_TYPE} ${WHICH_INTT} ${RUN_NUM}
mkdir -p $(dirname ${LIST})

ls -1 ${FILE} > ${LIST} 2>/dev/null
if [ ! -s ${LIST} ]; then
	echo -e "Expansion of ${FILE} failed\n"
	echo -e "Exiting\n"
	exit 1
	rm ${LIST}
fi

printf -v DATA ${DATA_FORMAT} ${RUN_NUM} ${WHICH_INTT}
mkdir -p $(dirname ${DATA})

# Macro
# root -q -b "${PWD}/macro.C(${RUN_NUM}, ${NUM_EVT}, \"${INTT_LIST}\")"
gdb -ex run --args root.exe "${PWD}/macro.C(${WHICH_INTT}, ${RUN_NUM}, ${NUM_EVT}, \"${LIST_FORMAT}\", \"${DATA_FORMAT}\")"

EXIT_VALUE="$?"

# Remove list files after finishing
rm -f ${LIST}

exit ${EXIT_VALUE}
