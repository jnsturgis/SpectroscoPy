# Documentation build.
#
#   make html     incremental -- fast, for writing
#   make clean    remove the build and the notebook execution cache
#   make rebuild  clean then html -- use before believing the navigation
#
# Why `rebuild` exists: Sphinx caches the table of contents between runs, so
# adding a page can leave the previous last page still claiming to be last, with
# no "next" link. An incremental build will not show it; a clean one will. CI
# always builds from scratch, so it never sees the stale version either.

SPHINXOPTS ?= -W --keep-going
SOURCEDIR  = .
BUILDDIR   = _build

.PHONY: html clean rebuild linkcheck

html:
	python3 -m sphinx -b html $(SPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/html

clean:
	rm -rf $(BUILDDIR) jupyter_execute

rebuild: clean html

linkcheck:
	python3 -m sphinx -b linkcheck $(SOURCEDIR) $(BUILDDIR)/linkcheck
