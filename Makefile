PYTHON ?= python3
MATH_DATABASE_DIR ?= ../math_database
DATA_DIR := $(CURDIR)/data
DOCS_DIR := $(CURDIR)/docs
LATEX_ROOT ?= $(HOME)/latex/PartialCubes
MIGRATION_DIR := $(CURDIR)/migration

.PHONY: validate build serve check clean migration-freeze migration-extract migration-init-queues migration-check

validate:
	$(PYTHON) scripts/validate_data.py --data-dir $(DATA_DIR)

build: validate
	$(PYTHON) $(MATH_DATABASE_DIR)/generate_website.py \
		--data_dir $(DATA_DIR) \
		--output_dir $(DOCS_DIR) \
		--deploy true
	$(PYTHON) scripts/normalize_generated_site.py --site-dir $(DOCS_DIR)
	touch $(DOCS_DIR)/.nojekyll

serve: validate
	$(PYTHON) $(MATH_DATABASE_DIR)/server.py --data_dir $(DATA_DIR)

check: build
	$(PYTHON) scripts/check_generated_site.py --site-dir $(DOCS_DIR)

migration-freeze:
	$(PYTHON) $(MIGRATION_DIR)/build_source_manifest.py --latex-root $(LATEX_ROOT)
	$(PYTHON) $(MIGRATION_DIR)/archive_dirty_sources.py --latex-root $(LATEX_ROOT)

migration-extract:
	$(PYTHON) $(MIGRATION_DIR)/extract_survey_atlas.py --latex-root $(LATEX_ROOT)
	$(PYTHON) $(MIGRATION_DIR)/extract_minor_inventory.py --latex-root $(LATEX_ROOT)
	$(PYTHON) $(MIGRATION_DIR)/index_latex_labels.py --latex-root $(LATEX_ROOT)
	$(PYTHON) $(MIGRATION_DIR)/audit_bibliography.py --latex-root $(LATEX_ROOT)
	$(PYTHON) $(MIGRATION_DIR)/audit_coverage.py --data-dir $(DATA_DIR)

migration-init-queues:
	$(PYTHON) $(MIGRATION_DIR)/initialize_review_queues.py

migration-check:
	$(PYTHON) $(MIGRATION_DIR)/extract_survey_atlas.py --latex-root $(LATEX_ROOT) --check
	$(PYTHON) $(MIGRATION_DIR)/extract_minor_inventory.py --latex-root $(LATEX_ROOT) --check
	$(PYTHON) $(MIGRATION_DIR)/index_latex_labels.py --latex-root $(LATEX_ROOT) --check
	$(PYTHON) $(MIGRATION_DIR)/check_source_manifest.py --latex-root $(LATEX_ROOT)
	$(PYTHON) $(MIGRATION_DIR)/audit_bibliography.py --latex-root $(LATEX_ROOT) --check
	$(PYTHON) $(MIGRATION_DIR)/audit_coverage.py --data-dir $(DATA_DIR) --check
	$(PYTHON) $(MIGRATION_DIR)/validate_review_queues.py

clean:
	rm -rf $(DOCS_DIR)
