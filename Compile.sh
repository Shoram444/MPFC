rm -rf build lib
rm FC.a

mkdir build lib

cd include

echo "rootcint: ./include/MPFeldman_Cousins.hh -> ./build/MPFeldman_Cousinsdict.cpp"
      rootcint -f ../build/MPFeldman_Cousinsdict.cpp   ./MPFeldman_Cousins.hh+

cd ../build

cmake ..
make

cd ..
cp FC.cpp ./build/
cd ./build/

root FC.cpp

##g++ -Wall -o FC.a -I./include -L./lib FC.cpp -lMPFC

##./FC.a
