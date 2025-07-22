include $(top_srcdir)/aminclude.am
ACLOCAL_AMFLAGS = -I m4

pkginclude_HEADERS = $(top_srcdir)/include/combinatorics.h \
$(top_srcdir)/include/permgrp.h $(top_srcdir)/include/perm.h \
$(top_srcdir)/include/permgrp_struct_p.h $(top_srcdir)/include/diophantine.h

lib_LTLIBRARIES = src/libcombinatorics.la
src_libcombinatorics_la_SOURCES = $(top_srcdir)/src/diophantine.c \
$(top_srcdir)/src/permgrp.c $(top_srcdir)/src/permgrp_aux_funcs.c \
$(top_srcdir)/src/permgrp_new.c
src_libcombinatorics_la_CFLAGS = -I$(top_srcdir)/include -Wall -Wextra -Wpedantic \
-Wno-unused-function
src_libcombinatorics_la_LIBADD = @LIBS@

if DOXY_
all-local:
	$(MAKE) doxygen-doc
endif

EXTRA_DIST = $(top_srcdir)/docs_pages
