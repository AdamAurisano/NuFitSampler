#!/usr/bin/env bash
# Download the NuFIT 6.0 chi^2 tables (IC19 w/o SK-atm = TBoff, IC24 with SK-atm = TByes)
# into data/nufit/ and verify them against the checksums recorded when this script was written.
# nu-fit.org currently serves an expired TLS certificate, hence the plain-http URL + checksum check.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dest="${here}/nufit"
mkdir -p "${dest}"
base="http://www.nu-fit.org/sites/default/files"

declare -A sha=(
  [v60.release-TBoff-IO.txt.xz]=1ce0c07ef6c4705d0e7ba27290cf011a56ad5f45a889fb73523b13b82000b713
  [v60.release-TBoff-NO.txt.xz]=575d96a0766b095a948780fbb2260c973d538f703ff05c65f7398a610013958e
  [v60.release-TByes-IO.txt.xz]=097cce0c1045863e38a302cb787c59162078f3ee881edd67172e0e4838cbfe02
  [v60.release-TByes-NO.txt.xz]=852b7791b276fe46749f00478bdd524426bae4fa42e9b5536687b001f2f11b14
)

for f in "${!sha[@]}"; do
  if [[ ! -f "${dest}/${f%.xz}" ]]; then
    curl -sSfL -o "${dest}/${f}" "${base}/${f}"
    echo "${sha[$f]}  ${dest}/${f}" | sha256sum -c --quiet -
    xz -d -f "${dest}/${f}"
  fi
  echo "ok ${dest}/${f%.xz}"
done
curl -sSfL -o "${dest}/v60.release-notes.pdf" "${base}/v60.release-notes.pdf" || true
