all: pdf
pdf: paper.pdf

.PHONY: clean 

paper.pdf: paper.md paper.bib *.jpg
	mkdir -p build
	TMPDIR="${PWD}/build" JOURNAL=joss singularity run --userns --bind=${PWD}:/data \
	    docker://openjournals/paperdraft:latest
	rm -rf build

arxiv: paper.md paper.bib *.jpg
	mkdir -p build
	mkdir -p submission
	TMPDIR="${PWD}/build" JOURNAL=joss singularity exec --userns --bind=${PWD}:/data \
	    docker://openjournals/paperdraft:latest sh -c \
	    'paperdraft paper.md --output submission/paper.tex; cp /usr/local/share/openjournals/joss/logo.png submission/'
	sed -i 's:/usr/local/share/openjournals/joss/logo.png:logo.png:' submission/paper.tex
	cp 02_operator_application.jpg submission
	tar -cf submission.tar submission/
	rm -r submission
	rm -r build

clean:
	rm -f paper.pdf paper.tex
	rm -f paper.aux paper.bcf paper.log paper.out paper.run.xml
	rm -rf build/
	rm -rf submission.tar submission/
