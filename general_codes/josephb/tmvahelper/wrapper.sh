#! /bin/bash

USR="$(id -u -n)"
PWD="$(pwd)"

# The location of this shell script
PWD=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

gdb -ex run --args root.exe -q -b macro/master.C
