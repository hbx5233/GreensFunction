# Makefile

SRC = ./src 
PROG = $(SRC)/greenmain.m 


figs: $(PROG)
	matlab -nodesktop -nosplash -nodisplay < $(PROG) 

clean:
	rm  Figs/*.pdf