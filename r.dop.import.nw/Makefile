MODULE_TOPDIR = ../..

PGM = r.dop.import.nw

include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Python.make

default: script $(PGM).md $(PGM).html

$(PGM).md: README.md
	cp README.md $(PGM).md

$(PGM).html: $(PGM).md
	echo "Creating extensions html file..."
	pandoc -f markdown -t html $(PGM).md -o $(PGM).html
	sed -i 's+<br />+\n<br>+g' $(PGM).html
