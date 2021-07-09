all: newton-cooling.png

newton-cooling.pdf: newton-cooling.tex
	pdflatex $<

newton-cooling.png: newton-cooling.pdf
	convert -density 300 $< -quality 90 $@

clean:
	rm -f *.aux
	rm -f *.log 

clean-all: clean
	rm -f *.pdf
	rm -f *.png

