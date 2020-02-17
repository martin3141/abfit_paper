all:
	latexmk main.tex
	latexmk supporting_info.tex
	
#	latexmk --pdf main.tex
#	latexmk --pdf supporting_info.tex

clean:
#   latexmk --pdf -C main.tex
	latexmk -C main.tex
	rm -f main.fff
	rm -f main.bbl
	rm -f main.synctex.gz
	rm -rf main.prv
	rm -rf auto
#   latexmk --pdf -C supporting_info.tex
	latexmk -C supporting_info.tex
	rm -f supporting_info.fff
	rm -f supporting_info.bbl
	rm -f supporting_info.synctex.gz
	rm -rf supporting_info.prv
