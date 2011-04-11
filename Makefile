RM	=	/bin/rm -f
RMDIR	=	/bin/rm -rf

all pdf ps:
	@echo The manual build system has been replaced by CMake. Please read the
	@echo README for instructions. You can save yourself some pain by running
	@echo 'make clean' now to get rid of old files that might break CMake
	@echo dependencies.

clean:
		$(RM) *.log *.lof *.lot *.bbl *.blg *.toc *.dvi *.aux *.ps *~ \#*\# *.idx *.ilg *.ind *.out *.pdf
		$(RM) progman.tex proglist.tex g_options.tex mdp_opt.tex programs.txt mdp_opt.html
		$(RMDIR) progman
