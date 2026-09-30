PYTHON ?= python3
MATH_DATABASE_DIR ?= ../math_database
DATA_DIR := $(CURDIR)/data
DOCS_DIR := $(CURDIR)/docs

.PHONY: validate build serve check clean

validate:
	$(PYTHON) scripts/validate_data.py --data-dir $(DATA_DIR)

build: validate
	$(PYTHON) $(MATH_DATABASE_DIR)/generate_website.py \
		--data_dir $(DATA_DIR) \
		--output_dir $(DOCS_DIR) \
		--deploy true
	touch $(DOCS_DIR)/.nojekyll

serve: validate
	$(PYTHON) $(MATH_DATABASE_DIR)/server.py --data_dir $(DATA_DIR)

check: build
	$(PYTHON) scripts/check_generated_site.py --site-dir $(DOCS_DIR)

clean:
	rm -rf $(DOCS_DIR)
