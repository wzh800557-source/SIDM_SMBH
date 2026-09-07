PYTHON ?= python3

.PHONY: test audit fluid-smoke gnc-prepare

test:
	MPLBACKEND=Agg $(PYTHON) tools/run_tests.py

audit:
	$(PYTHON) tools/audit_repository.py

fluid-smoke:
	PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg $(PYTHON) tools/smoke_test_fluid.py

gnc-prepare:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tools/build_gnc.py --prepare-only
