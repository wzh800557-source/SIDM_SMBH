#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from combine_ej_operators import combine  # noqa: E402


def write_operator(
    path: Path, seed: int, rates: np.ndarray, keep_second: bool,
    x_edges: np.ndarray | None = None,
) -> None:
    metadata = {
        "schema": "finite-angle-ej-jump-operator-v1",
        "profile_sha256": "profile",
        "df_table_sha256": "df",
        "mbh_msun": 4.0e6,
        "rh_pc": 20.0,
        "r_inner_pc": 1.0e-6,
        "reservoir_radius_pc": 2.0,
        "sigma0_over_m_cm2_g": 100.0,
        "w_kms": 80.0,
        "radial_bins": 8,
        "pairs_per_bin": 1000,
        "seed": seed,
        "energy_bins": 1,
        "angular_bins": 1,
        "domain_states_per_cell": 2,
        "state_count": 2,
        "post_capture_sentinel": -1,
        "post_reservoir_sentinel": -2,
        "reservoir_state_rule": "exact Kepler apocentre >= reporting radius",
        "nonkepler_exchange": "outer-fluid reservoir",
    }
    count = 2 if keep_second else 1
    np.savez_compressed(
        path,
        source1=np.array([1, 0], dtype=np.int64)[:count],
        source2=np.array([1, 1], dtype=np.int64)[:count],
        pre=np.array([1, 0], dtype=np.int64)[:count],
        post=np.array([0, -1], dtype=np.int64)[:count],
        rate_direct_msun_per_myr=rates[:count],
        rate_immediate_msun_per_myr=rates[:count],
        rate_direct_binding_msun_kms2_per_myr=(10.0 * rates[:count]),
        rate_immediate_binding_msun_kms2_per_myr=(10.0 * rates[:count]),
        sample_count=np.ones(count, dtype=np.int64),
        x_edges=(np.array([0.1, 1.0, np.inf]) if x_edges is None else x_edges),
        j_over_jlc_edges=np.array([1.0, 2.0, np.inf]),
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        first = directory / "first.npz"
        second = directory / "second.npz"
        write_operator(first, 11, np.array([2.0, 4.0]), True)
        write_operator(second, 22, np.array([4.0, 9.0]), False)
        arrays, metadata = combine([first, second])
        injection = (
            (arrays["pre"] == 1) & (arrays["post"] == 0)
        )
        capture = arrays["post"] == -1
        assert abs(float(np.sum(arrays["rate_direct_msun_per_myr"][injection])) - 3.0) < 1.0e-14
        # The second realization did not sample the capture record. Its zero
        # contribution must remain in the denominator, giving (4+0)/2=2.
        assert abs(float(np.sum(arrays["rate_direct_msun_per_myr"][capture])) - 2.0) < 1.0e-14
        assert abs(float(np.sum(
            arrays["rate_direct_binding_msun_kms2_per_myr"][capture]
        )) - 20.0) < 1.0e-14
        assert metadata["capture_energy_columns_present"] is True
        assert metadata["seeds"] == [11, 22]
        assert metadata["input_paths"] == ["first.npz", "second.npz"]
        assert metadata["input_edge_grid_max_relative_spread"]["x_edges"] == 0.0

        # Transcendental grid construction can differ by an ulp across CPU
        # types. That roundoff is accepted and recorded, while a physical grid
        # change remains an error.
        roundoff = directory / "roundoff.npz"
        write_operator(
            roundoff, 33, np.array([4.0, 9.0]), False,
            np.array([np.nextafter(0.1, 1.0), 1.0, np.inf]),
        )
        _, roundoff_metadata = combine([first, roundoff])
        assert 0.0 < roundoff_metadata[
            "input_edge_grid_max_relative_spread"
        ]["x_edges"] < 1.0e-13
        changed = directory / "changed.npz"
        write_operator(
            changed, 44, np.array([4.0, 9.0]), False,
            np.array([0.1001, 1.0, np.inf]),
        )
        try:
            combine([first, changed])
        except RuntimeError as error:
            assert "physically different x_edges" in str(error)
        else:
            raise AssertionError("the operator combiner accepted a changed grid")
        changed_sentinel = directory / "changed_sentinel.npz"
        write_operator(
            changed_sentinel, 55, np.array([4.0, 9.0]), False,
            np.array([0.1, 1.0, 3.0]),
        )
        try:
            combine([first, changed_sentinel])
        except RuntimeError as error:
            assert "finite/infinite x_edges" in str(error)
        else:
            raise AssertionError("the operator combiner accepted a changed sentinel")
        duplicate = directory / "duplicate.npz"
        write_operator(duplicate, 11, np.array([4.0, 9.0]), False)
        try:
            combine([first, duplicate])
        except RuntimeError as error:
            assert "repeated EJ random seeds" in str(error)
        else:
            raise AssertionError("the operator combiner accepted a repeated seed")
    print("PASS: zero-filled combination of independent EJ operators")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
