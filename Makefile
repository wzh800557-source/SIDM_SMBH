PYTHON ?= python3

.PHONY: test audit

test:
	MPLBACKEND=Agg $(PYTHON) tools/run_tests.py

audit:
	$(PYTHON) tools/audit_repository.py
