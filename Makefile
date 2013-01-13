package = pear
version = 0.0.1
tarname = $(package)
distdir = $(tarname)-$(version)

all clean install pear:
	$(MAKE) -C src $@

dist: $(distdir).tar.gz

$(distdir).tar.gz: FORCE $(distdir)
	tar chof - $(distdir) |\
	  gzip -9 -c >$(distdir).tar.gz
	rm -rf $(distdir)

$(distdir):
	mkdir -p $(distdir)/src
	cp Makefile $(distdir)
	cp README $(distdir)
	cp src/Makefile $(distdir)/src
	cp src/pear.c $(distdir)/src
	cp src/fastq.c $(distdir)/src
	cp src/fastq.h $(distdir)/src

distcheck: $(distdir).tar.gz
	gzip -cd $+ | tar xvf -
	$(MAKE) -C $(distdir) all clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz\
	  ready for distribution."

FORCE:
	-rm $(distdir).tar.gz &> /dev/null
	-rm -rf $(distdir) &> /dev/null

.PHONY: FORCE all clean dist distcheck
