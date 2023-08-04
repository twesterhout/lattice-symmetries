#!/bin/bash

set -e
set -o pipefail

show_help() {
    cat <<-EOF
Usage: ${0##*/} [-h] [NAME]

    -h          display this help and exit
EOF
}

OPTIND=1
while getopts h: opt; do
    case $opt in
    h)
        show_help
        exit 0
        ;;
    *)
        show_help >&2
        exit 1
        ;;
    esac
done
shift "$((OPTIND - 1))"

PROJECT_NAME="lattice-symmetries-kernel"
if [[ $# -ge 1 ]]; then
    PROJECT_NAME=$1
fi

# NOTE: Determine the location of the script.
# See https://stackoverflow.com/a/246128 for details
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

cd "$SCRIPT_DIR"
# Build the container
CONTAINER_PATH=$(nix build --no-link --print-out-paths .)
# Copy the container to the current directory
nix shell nixpkgs#coreutils -c sh -c "cp -b -f -v '$CONTAINER_PATH' container.img"
# Install the ipython kernel
singularity exec container.img python3 -m ipykernel install \
    --user \
    --name="$PROJECT_NAME"

cat >"$HOME/.local/share/jupyter/kernels/$PROJECT_NAME/kernel.json" <<EOF
{
 "argv": [
  "singularity",
  "exec",
  "$(realpath ./container.img)",
  "python3",
  "-m",
  "ipykernel_launcher",
  "-f",
  "{connection_file}"
 ],
 "display_name": "$PROJECT_NAME",
 "language": "python",
 "metadata": {
  "debugger": false
 }
}
EOF
