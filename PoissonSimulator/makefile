simulator: Main.o alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o fasttransforms.o integration.o interpolation.o linalg.o optimization.o solvers.o specialfunctions.o statistics.o stdafx.o Poisson1d.o Poisson2d.o
	g++ -o simulator Main.o alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o fasttransforms.o integration.o interpolation.o linalg.o optimization.o solvers.o specialfunctions.o statistics.o stdafx.o Poisson1d.o Poisson2d.o
     
Main.o: Main.cpp 
	g++ -c $(CPPFLAGS) Main.cpp

alglibinternal.o: alglibinternal.cpp alglibinternal.h
	g++ -c $(CPPFLAGS) alglibinternal.cpp
	
alglibmisc.o:	alglibmisc.cpp alglibmisc.h
	g++ -c $(CPPFLAGS) alglibmisc.cpp
	
ap.o: ap.cpp ap.h
	g++ -c $(CPPFLAGS) ap.cpp

dataanalysis.o: dataanalysis.cpp dataanalysis.h
	g++ -c $(CPPFLAGS) dataanalysis.cpp
	
diffequations.o: diffequations.cpp diffequations.h
	g++ -c $(CPPFLAGS) diffequations.cpp
	
fasttransforms.o: fasttransforms.cpp fasttransforms.h
	g++ -c $(CPPFLAGS) fasttransforms.cpp
	
integration.o: integration.cpp integration.h
	g++ -c $(CPPFLAGS) integration.cpp
	
interpolation.o: interpolation.cpp interpolation.h
	g++ -c $(CPPFLAGS) interpolation.cpp
	
linalg.o: linalg.cpp linalg.h
	g++ -c $(CPPFLAGS) linalg.cpp
	
optimization.o: optimization.cpp optimization.h
	g++ -c $(CPPFLAGS) optimization.cpp
	
solvers.o: solvers.cpp solvers.h
	g++ -c $(CPPFLAGS) solvers.cpp
	
specialfunctions.o: specialfunctions.cpp specialfunctions.h
	g++ -c $(CPPFLAGS) specialfunctions.cpp
	
statistics.o: statistics.cpp statistics.h
	g++ -c $(CPPFLAGS) statistics.cpp
	
stdafx.o: stdafx.cpp stdafx.h
	g++ -c $(CPPFLAGS) stdafx.cpp
	
Poisson1d.o: Poisson1d.cpp Poisson1d.hpp
	g++ -c $(CPPFLAGS) Poisson1d.cpp
	
Poisson2d.o: Poisson2d.cpp Poisson2d.hpp
	g++ -c $(CPPFLAGS) Poisson2d.cpp
	

clean:
	rm simulator Main.o alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o fasttransforms.o integration.o interpolation.o linalg.o optimization.o solvers.o specialfunctions.o statistics.o stdafx.o Poisson1d.o Poisson2d.o