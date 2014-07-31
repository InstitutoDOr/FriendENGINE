include $(FSLCONFDIR)/default.mk

PROJNAME = fslpres

LIBS = -lss_32R -lm

SCRIPTS = fslpres eqn2gif makehtmeqns

all:	imagesinstall

imagesinstall:
	/bin/rm -rf $(DESTDIR)/etc/fslpres ; /bin/cp -r lib $(DESTDIR)/etc/fslpres;
	/bin/cp .latex2html-init $(DESTDIR)/etc/fslpres

