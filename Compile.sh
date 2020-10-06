	rm -rf build lib ./src/dicts/
	mkdir build lib

echo "#############################################"
echo "#         GENERATE ROOT DICTIONARIES        #"
echo "#############################################"
echo " "

	cd include

	echo "rootcint: ./include/MPFeldman_Cousins.hh -> ./src/dicts/MPFeldman_Cousinsdict.cpp + ./lib/MPFeldman_Cousinsdict_rdict.pcm"
      	rootcint -f ../lib/MPFeldman_Cousinsdict.cpp   MPFeldman_Cousins.hh+
	
	echo " "
	echo "Dictionaries generated!"
	echo " "

	cd ..

	cp -R ./lib/ ./src/dicts/
	rm -rf ./lib/*.cpp
	rm -rf ./src/dicts/*.pcm

echo "#############################################"
echo "#         MODULE COMPILATION (CMAKE)        #"
echo "#############################################"
echo " "	

	cd build

	cmake .. || { echo $'\n****Cmake failed, installation of MPFC library aborted!****' ; exit 1; }

echo " "
echo "#############################################"
echo "#         MODULE COMPILATION (MAKE)         #"
echo "#############################################"
echo " "

	make

	cp libMPFC.so ../lib/libMPFC.so
	rm libMPFC.so

