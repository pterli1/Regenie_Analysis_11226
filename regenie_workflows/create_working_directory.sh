#!/bin/bash

set -euo pipefail

# Helper script to set up fedcat from cloned repo
# Usage: create_working_directory.sh [target_directory (default: fedcat_runs)] [repo_root (default: this script's repo)]

target_dir=${1:-fedcat_runs}
repo_root=${2:-"$(dirname "${BASH_SOURCE[0]}")/.."}

# If you run this script while another venv is active, `pip` may point to that venv
# and dependencies may end up installed there (not desired).
if [[ -n "${VIRTUAL_ENV:-}" ]]; then
    echo "Error: Detected active virtualenv: ${VIRTUAL_ENV}"
    echo "  For dependencies to be installed into the working directory venv, run this script from a shell with no venv active."
    echo "  (In your current shell, run: deactivate)"
    exit 1
fi

# Normalise repo_root for safer path joins
repo_root="$(cd "$repo_root" && pwd)"
fedcat_dir="$repo_root/regenie_workflows/federated_analysis"
fedcat_rs_src="$fedcat_dir/fedcat_rs"

# fedcat's fedcat_rs requirement, e.g. "fedcat_rs>=0.1.5"
required_spec="$(sed -n 's/^[[:space:]]*"\(fedcat_rs[^"]*\)".*/\1/p' "$fedcat_dir/pyproject.toml" | head -n1)"
if [[ -z "$required_spec" ]]; then
    echo "Error: could not find a fedcat_rs requirement in $fedcat_dir/pyproject.toml" >&2
    exit 1
fi

# Download the newest fedcat_rs wheel from the public GHCR package (pushed by the
# fedcat-meta-docker action, Linux x86_64 only) into $1
fetch_ghcr_wheel() {
    if [[ "$(uname -s)/$(uname -m)" != "Linux/x86_64" ]]; then
        echo "WARNING: prebuilt fedcat_rs wheels are not yet available for $(uname -s)/$(uname -m) (Linux x86_64 only)" >&2
        return 1
    fi
    python3 - "$1" <<'PY'
import hashlib
import json
import pathlib
import re
import sys
import urllib.request

name = "exetergenetics/biobank_workflows/fedcat_rs"
out_dir = pathlib.Path(sys.argv[1])


def get(path, accept="application/json"):
    headers = {"Authorization": f"Bearer {token}", "Accept": accept}
    with urllib.request.urlopen(urllib.request.Request(f"https://ghcr.io/v2/{name}/{path}", headers=headers)) as resp:
        return resp.read()


try:
    with urllib.request.urlopen(f"https://ghcr.io/token?scope=repository:{name}:pull") as resp:
        token = json.load(resp)["token"]
    # Tags are fedcat_rs versions, each holding a single wheel
    tags = [t for t in json.loads(get("tags/list"))["tags"] if re.fullmatch(r"[0-9.]+", t)]
    latest = max(tags, key=lambda t: tuple(map(int, t.split("."))))
    layer = json.loads(get(f"manifests/{latest}", "application/vnd.oci.image.manifest.v1+json"))["layers"][0]
    wheel = get(f"blobs/{layer['digest']}")
    if "sha256:" + hashlib.sha256(wheel).hexdigest() != layer["digest"]:
        raise ValueError("checksum mismatch")
except Exception as e:
    sys.exit(f"Could not fetch fedcat_rs wheel from ghcr.io/{name}: {e}")

filename = layer["annotations"]["org.opencontainers.image.title"]
(out_dir / filename).write_bytes(wheel)
print(f"Downloaded {filename} from ghcr.io/{name}")
PY
}

mkdir -p "$target_dir"
cd "$target_dir"

cp "$fedcat_dir/setup_files/variant_density_blocks.tsv" "./"

export PIP_USER=false
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip wheel

# fedcat_rs is not on PyPI, so install it before fedcat. Try the GHCR wheel and any
# local wheels from build_local_wheel.sh, falling back to a slow source build.
wheel_dir="$(mktemp -d)"
trap 'rm -rf "$wheel_dir"' EXIT
fetch_ghcr_wheel "$wheel_dir" || true
if ! pip install --no-index --find-links "$wheel_dir" --find-links "$fedcat_rs_src/dist" "$required_spec"; then
    echo "WARNING: no prebuilt fedcat_rs wheel satisfies $required_spec, building from source instead." >&2
    echo "  This needs a Rust toolchain (https://rustup.rs) and can take several minutes." >&2
    pip install "$fedcat_rs_src"
fi

pip install -e "$fedcat_dir"

echo "Working directory setup complete at $target_dir"
echo "Activate the virtual environment with: source $target_dir/venv/bin/activate"
echo "Then run: fedcat discover --config config.json"
