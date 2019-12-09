all:
	latexmk --pdf main.tex

clean:
	latexmk --pdf -C main.tex
	rm -f main.fff
	rm -f main.bbl
	rm -f main.synctex.gz
	rm -rf auto
	rm -rf main.prv
