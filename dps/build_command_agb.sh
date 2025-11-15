#!/usr/bin/env -S bash --login
set -euo pipefail
basedir=$( cd "$(dirname "$0")" ; pwd -P )
conda env update -f ${basedir}/above_env_r_3.1.4.yml





