#Makefile

espic:define.hpp ConstParam_pulse.cpp main.cpp SolveField.cpp SetParticle.cpp SolveMotion.cpp OutPut.cpp

	mpic++ -std=c++11 -o Pulse ConstParam_pulse.cpp main.cpp SolveField.cpp SetParticle.cpp  SolveMotion.cpp OutPut.cpp
	
	$(shell mkdir -p Electron_Data) $(shell mkdir -p Ion_Data) $(shell mkdir -p Grid_Data)

clean:
	rm -f *.gif *.dat
	rm -f Electron_Data/* Ion_Data/* Grid_Data/*

cleanall:
	rm -r Pulse*
	rm -f *.gif *.dat
	rm -r Electron_Data Ion_Data Grid_Data
	
