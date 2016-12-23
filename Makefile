# Makefile
 
PROG = ./src/greenmain.m 

figs: $(PROG)
	mkdir Figs
	matlab -nodesktop -nosplash -nodisplay < $(PROG) 

clean:
	rm  -rf Figs