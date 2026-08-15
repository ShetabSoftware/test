#!/usr/bin/env bash
# =====================================================================
#  run_sim.sh  -  co-simulate the RTL against the MATLAB golden model.
#
#  Regenerates the golden vectors, analyses the whole design, then runs
#  every stage testbench.  A stage passes only if its output is
#  BIT-IDENTICAL to the model; there is no tolerance anywhere.
#
#    ./rtl/sim/run_sim.sh            # all stages, 2 ms of data
#    ./rtl/sim/run_sim.sh ddc fir    # only the named stages
#    MS=6 ./rtl/sim/run_sim.sh       # longer run
#
#  Requires: ghdl (>=3.0) and octave-cli or matlab.
# =====================================================================
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="$ROOT/build"
GOLD="$BUILD/gold"
MS="${MS:-2}"
STD=93

mkdir -p "$BUILD"
cd "$BUILD" || exit 1

# ---------------------------------------------------------------- vectors
if [ ! -f "$GOLD/s01_ddc_out.txt" ] || [ "${REGEN:-0}" = "1" ]; then
  echo "=== generating golden vectors (${MS} ms) ==="
  ( cd "$ROOT/matlab" && octave-cli --no-gui -q --eval \
      "asp_startup; asp_export_vhdl; asp_golden_model('durationMs',${MS},'outDir','${GOLD}','verbose',false);" ) \
    || { echo "golden model failed"; exit 1; }
fi

# ---------------------------------------------------------------- analyse
echo "=== analysing RTL ==="
SRC=(
  "$ROOT/rtl/pkg/asp_pkg.vhd"
  "$ROOT/rtl/pkg/asp_coef_pkg.vhd"
  "$ROOT/rtl/core/asp_ddc_mixer.vhd"
  "$ROOT/rtl/core/asp_hb_decim2.vhd"
  "$ROOT/rtl/core/asp_fir_shape.vhd"
  "$ROOT/rtl/core/asp_cov_accum.vhd"
  "$ROOT/rtl/core/asp_rsqrt.vhd"
  "$ROOT/rtl/core/asp_isqrt.vhd"
  "$ROOT/rtl/core/asp_whiten.vhd"
  "$ROOT/rtl/core/asp_cordic.vhd"
  "$ROOT/rtl/core/asp_jacobi_evd.vhd"
  "$ROOT/rtl/core/asp_detect.vhd"
  "$ROOT/rtl/core/asp_weight_calc.vhd"
  "$ROOT/rtl/core/asp_beamformer.vhd"
  "$ROOT/rtl/core/asp_tx_scale.vhd"
  "$ROOT/rtl/top/asp_dwell_engine.vhd"
  "$ROOT/rtl/top/asp_datapath.vhd"
  "$ROOT/rtl/top/asp_axi_lite_regs.vhd"
  "$ROOT/rtl/top/asp_top.vhd"
  "$ROOT/rtl/tb/asp_tb_pkg.vhd"
)
for f in "${SRC[@]}"; do
  [ -f "$f" ] || continue
  ghdl -a --workdir=. --std=$STD "$f" || { echo "ANALYSE FAILED: $f"; exit 1; }
done
for f in "$ROOT"/rtl/tb/tb_*.vhd; do
  [ -f "$f" ] || continue
  ghdl -a --workdir=. --std=$STD "$f" || { echo "ANALYSE FAILED: $f"; exit 1; }
done

# ---------------------------------------------------------------- run
declare -A TBS=(
  [ddc]="tb_asp_ddc_mixer"
  [hb]="tb_asp_hb_decim2"
  [fir]="tb_asp_fir_shape"
  [cov]="tb_asp_cov_accum"
  [rsqrt]="tb_asp_rsqrt"
  [whiten]="tb_asp_whiten"
  [cordic]="tb_asp_cordic"
  [evd]="tb_asp_jacobi_evd"
  [detect]="tb_asp_detect"
  [weights]="tb_asp_weight_calc"
  [beam]="tb_asp_beamformer"
  [tx]="tb_asp_tx_scale"
  [chain]="tb_asp_datapath"
)
# Fast blocks first so a regression is reported in seconds rather than
# after the multi-minute streaming tests.
ORDER=(rsqrt cordic detect whiten evd weights cov ddc hb tx fir beam chain)

WANT=("$@")
[ ${#WANT[@]} -eq 0 ] && WANT=("${ORDER[@]}")

fail=0
run=0
for key in "${ORDER[@]}"; do
  for w in "${WANT[@]}"; do
    if [ "$key" = "$w" ]; then
      tb="${TBS[$key]}"
      if ghdl -e --workdir=. --std=$STD "$tb" >/dev/null 2>&1; then
        echo "=== $tb ==="
        run=$((run+1))
        if ! ghdl -r --workdir=. --std=$STD "$tb" 2>&1 | grep -Ev '^$'; then
          fail=$((fail+1))
        fi
      else
        echo "--- $tb not present, skipped"
      fi
    fi
  done
done

echo
if [ $fail -eq 0 ]; then
  echo "==================================================="
  echo " ALL $run TESTBENCHES PASSED (bit-exact vs MATLAB)"
  echo "==================================================="
else
  echo "*** $fail of $run TESTBENCHES FAILED ***"
fi
exit $fail
