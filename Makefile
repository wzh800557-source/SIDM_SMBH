PYTHON ?= python3

.PHONY: test audit fluid-smoke gnc-prepare calibrate final-calibration

test:
	MPLBACKEND=Agg $(PYTHON) tools/run_tests.py

audit:
	$(PYTHON) tools/audit_repository.py

fluid-smoke:
	PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg $(PYTHON) tools/smoke_test_fluid.py

gnc-prepare:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tools/build_gnc.py --prepare-only

calibrate:
	PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg $(PYTHON) fp_solver/calibrate_cross_regime.py --out validation/cross_regime_calibration.json
	PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg $(PYTHON) figures/plot_cross_regime_calibration.py --outdir validation/cross_regime_figures

final-calibration: calibrate
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tools/build_final_calibration_audit.py --out validation/final_calibration_audit.json
	PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg $(PYTHON) tools/run_tests.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tools/audit_repository.py
