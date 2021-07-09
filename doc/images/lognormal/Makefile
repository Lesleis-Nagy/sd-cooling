all: lognormal.png

lognormal.pdf: lognormal.tex
	pdflatex $<

lognormal.png: lognormal.pdf
	convert -density 300 $< -quality 90 $@

clean:
	rm -f *.aux
	rm -f *.log 

clean-all: clean
	rm -f *.pdf
	rm -f *.png

