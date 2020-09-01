rm FC.a

g++ -c -Wall -o lib/MPFeldman_Cousins.o -I./include ./src/MPFeldman_Cousins.cpp 
ar rvs lib/libMPFeldman_Cousins.a lib/MPFeldman_Cousins.o

g++ -Wall -o FC.a -I./include -L./lib FC.cpp -lMPFeldman_Cousins

./FC.a