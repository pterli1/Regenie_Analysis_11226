#!/bin/bash

set -euo pipefail

# Build a local fedcat_rs wheel for faster private installs.
# Output wheel(s) are written to fedcat_rs/dist/.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

python3 -m pip install --upgrade pip wheel maturin
python3 -m pip wheel --wheel-dir "$script_dir/dist" "$script_dir"

echo "Built wheel(s) in: $script_dir/dist"
echo "Use them during setup via:"
echo "  pip install --find-links $script_dir/dist fedcat_rs==0.1.1"
