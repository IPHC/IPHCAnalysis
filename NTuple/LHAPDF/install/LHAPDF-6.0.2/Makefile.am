ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src include wrappers
DIST_SUBDIRS = $(SUBDIRS) examples tests doc

bin_SCRIPTS = lhapdf-config
dist_pkgdata_DATA = lhapdf.conf pdfsets.index


tests:
	(cd tests && $(MAKE))

examples:
	(cd examples && $(MAKE))

if WITH_DOXYGEN
doxy:
	(cd doc && $(MAKE) doxy)
endif

if WITH_PDFLATEX
manual:
	(cd doc && $(MAKE) manual)
endif


## Remove temporary program generated by AC compiler tests
mostlyclean-local:
	rm -rf a.out

## Clean out SVN and hg files
dist-hook:
	rm -rf `find $(distdir) -name ".svn" -o -name ".hg"`

## Remove the unused & misleading share/LHAPDF/sets directory on both install and uninstall!
# TODO: Remove this hook entirely after a couple of versions
install-data-hook:
	test -e $(pkgdatadir)/sets && rmdir $(pkgdatadir)/sets || true
uninstall-hook:
	test -e $(pkgdatadir)/sets && rmdir $(pkgdatadir)/sets || true
