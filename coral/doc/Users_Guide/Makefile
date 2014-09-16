#-------------------------------------
#   Main targets
#-------------------------------------
all       : users_guide.ps users_guide.pdf
install   : 
uninstall :
dist      :
clean     :
	rm -f core core.* *.o *~ *.log *.ii *.ti *.out *.dvi *.aux 
realclean : clean
	cd ..; rm -f *.pdf *.ps
distclean : clean

#-------------------------------------
#  suffix rules
#-------------------------------------
.SUFFIXES : .tex .dvi .ps .pdf .aux .epsi .ps2 .eps

.tex.aux :
	latex $*
.tex.dvi :
	latex $* ; latex $*
.dvi.ps :
	dvips -t letter -o $@ $<
.ps.pdf :
	ps2pdf $< $@
.ps.epsi :
	ps2ps $< temp.ps ; ps2epsi temp.ps $@ ; rm temp.ps
.ps.ps2 :
	psnup -2 -pletter -Pletter $< $@
.ps.eps :
	ps2ps $< $@
  
#-------------------------------------
#  specific dependencies
#-------------------------------------
users_guide.dvi : users_guide.tex globalparameters.tex inputfiles.tex \
    objects.tex actions.tex kernels.tex
