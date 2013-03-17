ifndef PROCS
	PROCS=8
endif

TRACEDIR= ./projections
#PROJ = 1
ifdef PROJ
	TRACEMODE=-tracemode projections
	TRACEROOT=+traceroot $(TRACEDIR)
endif

LIBDIR = ../../lib
SRCDIR = .
CHARM_ON = 1
COPTS = -g -DCHARM=$(CHARM_ON) -O2 #-ansi -pedantic -Wall

CHARMC = charmc $(OPTS) $(TRACEMODE)
CC = g++

INCLUDES = -I$(SRCDIR) -I$(LIBDIR)
ALL = tracker testtracker #trackerserial
VPATH = .:$(LIBDIR)

RUN = charmrun

.PHONY: all clean clobber

all: $(ALL)

clean:
	-rm -rf *.o *~ core

clobber:
	-rm -rf charmrun \#*\# *.o *~ *.decl.h *.def.h core core.* $(TRACEDIR) $(ALL)

tracker: Tracker.o TrackParticles.o Search.o Cost.o Input.o Regression.o fastRegression.o mergeRegression.o FrameSet.o OutputG.o Merger.o kinetics.o
	$(CHARMC) $(COPTS) $^ -o $@ -lm

testtracker: TestTracker.o Input.o Output.o Goldtest.o AnalyzeInput.o Cost.o Regression.o fastRegression.o mergeRegression.o
	$(CC) $(COPTS) $^ -o $@ -lm

#trackerserial: TrackerSerial.o Input.o Regression.o mergeRegression.o fastRegression.o Cost.o Search.o Output.o Goldtest.o linuxtimer.o TrackParticles.o AnalyzeInput.o
#	$(CC) $(COPTS) $(INCLUDES) $^ -o $@ -lm 
	
Tracker.decl.h: Main.ci
	$(CHARMC) $(COPTS) Main.ci

.cpp.o:
	$(CHARMC) $(COPTS) $(INCLUDES) -c $<

Tracker.o: Tracker.cpp Main.ci Tracker.decl.h interface.h FrameSet.h FrameSet.cpp OutputG.h OutputG.cpp Merger.cpp Merger.h Tracker.decl.h
OutputG.o: OutputG.h OutputG.cpp interface.h Tracker.decl.h
Merger.o: interface.h Merger.h Merger.cpp Cost.h Cost.cpp Tracker.decl.h
FrameSet.o: FrameSet.h FrameSet.cpp Tracker.cpp interface.h Tracker.decl.h
TrackParticles.o: TrackParticles.h TrackParticles.cpp interface.h Tracker.decl.h
Regression.o: Regression.h Regression.cpp interface.h 
Search.o: Search.h Search.cpp interface.h Tracker.decl.h
Cost.o: Cost.h Cost.cpp interface.h Tracker.decl.h
Input.o: Input.h Input.cpp interface.h Tracker.decl.h
fastRegression.o: fastRegression.cpp fastRegression.h interface.h Tracker.decl.h
mergeRegression.o: mergeRegression.cpp mergeRegression.h interface.h Tracker.decl.h
Goldtest.o: Goldtest.h Goldtest.cpp interface.h
TestTracker.o: TestTracker.cpp interface.h
TrackerSerial.o: TrackerSerial.cpp interface.h
linuxtimer.o: linuxtimer.cpp linuxtimer.h
AnalyzeInput.o: AnalyzeInput.h AnalyzeInput.cpp interface.h
kinetics.o: kinetics.h kinetics.cpp interface.h

run: tracker
ifdef PROJ
	-mkdir tracker-${PROCS}
endif
	-mkdir ../output
	$(RUN) +p${PROCS} ./$^ scaling/256p512f out -$(PROCS) $(TRACEROOT)
