.PHONY: zip clean

zip:
	rm -f xsladk07.zip
	mkdir zip
	cp calculators/BatchMandelCalculator.cc zip/BatchMandelCalculator.cc
	cp calculators/BatchMandelCalculator.h zip/BatchMandelCalculator.h
	cp calculators/LineMandelCalculator.cc zip/LineMandelCalculator.cc
	cp calculators/LineMandelCalculator.h zip/LineMandelCalculator.h
	cp eval.png zip/eval.png
	cp MB-xsladk07.txt zip/MB-xsladk07.txt
	cd zip && zip xsladk07.zip BatchMandelCalculator.cc BatchMandelCalculator.h LineMandelCalculator.cc LineMandelCalculator.h eval.png MB-xsladk07.txt
	mv zip/xsladk07.zip xsladk07.zip
	rm -r zip