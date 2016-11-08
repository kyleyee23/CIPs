#
# target entry to build program executable from program and mylib 
# object files 
#
program: CIP_lensing.cpp wignerSymbols/src/wignerSymbols-cpp.cpp
	g++ -o test_CIP_lensing CIP_lensing.cpp globaldata.h wignerSymbols/src/wignerSymbols-cpp.cpp
