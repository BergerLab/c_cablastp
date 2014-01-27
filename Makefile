all: cablastp-compress cablastp-decompress cablastp-search

cablastp-compress: src/cablastp-compress
	cp src/cablastp-compress .

src/cablastp-compress:
	(cd src && make compress)

cablastp-decompress: src/cablastp-decompress
	cp src/cablastp-decompress .

src/cablastp-decompress:
	(cd src && make decompress)

cablastp-search: src/cablastp-search
	cp src/cablastp-search .

src/cablastp-search:
	(cd src && make search)


clean:
	(cd src && make clean)
	rm -f cablastp-*

push:
	git push origin master
	git push github master

