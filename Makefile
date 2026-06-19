# Makefile — molecular models for the "Ecosystems as energy fields" figures.
# All Python targets assume the `energetics-pymol` conda env (see `make env`).

ENV    := energetics-pymol
RUN    := conda run -n $(ENV)
PY     := $(RUN) python

.PHONY: all env pdb render thumbs catalog lock clean distclean help

help:
	@echo "make env       create the conda env (rdkit + openbabel + pymol-open-source)"
	@echo "make all       build every .pdb, render every .png, write the contact sheet"
	@echo "make pdb       build pdb/*.pdb from the manifest (RDKit)"
	@echo "make render    write scripts/*.pml and ray-trace renders/*.png (PyMOL)"
	@echo "make thumbs    fast low-res verification renders"
	@echo "make catalog   print the molecule catalog (no RDKit needed)"
	@echo "make lock      export environment.yml from the live env"
	@echo "make clean     remove generated pdb/, scripts/, renders/"

env:
	conda create -n $(ENV) -c conda-forge --strict-channel-priority -y \
		python=3.11 rdkit openbabel pymol-open-source numpy pyyaml pillow
	@echo "Activate with:  conda activate $(ENV)"

all:
	$(PY) build_all.py

pdb:
	$(PY) build_pdb.py

render:
	$(PY) render.py

thumbs:
	$(PY) render.py --thumb

catalog:
	$(PY) molecules.py

lock:
	conda env export -n $(ENV) --no-builds > environment.yml
	@echo "wrote environment.yml"

clean:
	rm -rf pdb scripts renders

distclean: clean
	conda env remove -n $(ENV) -y || true
