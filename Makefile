all: figures/fig1.eps figures/fig2.eps figures/fig3.eps figures/fig4.eps figures/fig5.eps figures/fig6.eps figures/fig7.eps figures/fig8.eps figures/fig9.eps figures/fig10.eps 
	latexmk -pdf main.tex
	latexmk -pdf supporting_info.tex

figures/fig1.eps: figures/fig1.tiff
	convert figures/fig1.tiff eps3:figures/fig1.eps

figures/fig2.eps: figures/fig2.tiff
	convert figures/fig2.tiff eps3:figures/fig2.eps

figures/fig3.eps: figures/fig3.tiff
	convert figures/fig3.tiff eps3:figures/fig3.eps

figures/fig4.eps: figures/fig4.tiff
	convert figures/fig4.tiff eps3:figures/fig4.eps

figures/fig5.eps: figures/fig5.tiff
	convert figures/fig5.tiff eps3:figures/fig5.eps

figures/fig6.eps: figures/fig6.tiff
	convert figures/fig6.tiff eps3:figures/fig6.eps

figures/fig7.eps: figures/fig7.tiff
	convert figures/fig7.tiff eps3:figures/fig7.eps

figures/fig8.eps: figures/fig8.tiff
	convert figures/fig8.tiff eps3:figures/fig8.eps

figures/fig9.eps: figures/fig9.tiff
	convert figures/fig9.tiff eps3:figures/fig9.eps

figures/fig10.eps: figures/fig10.tiff
	convert figures/fig10.tiff eps3:figures/fig10.eps

clean:
	latexmk -pdf -C main.tex
	rm -f main.fff
	rm -f main.bbl
	rm -f main.synctex.gz
	rm -rf main.prv
	rm -rf auto

	latexmk -pdf -C supporting_info.tex
	rm -f supporting_info.fff
	rm -f supporting_info.bbl
	rm -f supporting_info.synctex.gz
	rm -rf supporting_info.prv
