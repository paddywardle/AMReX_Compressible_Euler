myProgram: main.o EulerSystem.o EulerSystem.h NumericalMethod.o NumericalMethod.h SettingsData.o InitialCondTests.o
	g++ -I/usr/include/eigen3/ main.o EulerSystem.o NumericalMethod.o SettingsData.o InitialCondTests.o -lconfig++ -o myProgram -O3

main.o: main.cpp
	g++ -c main.cpp -o main.o -O3

EulerSystem.o: EulerSystem.cpp EulerSystem.h
	g++ -c EulerSystem.cpp -o EulerSystem.o -O3

NumericalMethod.o: NumericalMethod.cpp NumericalMethod.h
	g++ -c NumericalMethod.cpp -o NumericalMethod.o -O3

InitialCondTests.o: InitialCondTests.cpp InitialCondTests.h
	g++ -c InitialCondTests.cpp -o InitialCondTest.o -O3

SettingsData.o: SettingsData.cpp SettingsData.h
	g++ -c SettingsData.cpp -o SettingsData.o -O3

