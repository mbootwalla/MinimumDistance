## A Makefile for running the cbs vignette
## must run make from a fat node
Beaty : Beaty.Rnw
	echo "Stangle(\"Beaty.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/Beaty.R/" > Beaty.R.sh
	qsub -m e -r y -cwd -l mem_free=5G,h_vmem=12G Beaty.R.sh

segmentation: segmentation.Rnw
	echo "Stangle(\"segmentation.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/segmentation.R/" > segmentation.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=5G,h_vmem=12G segmentation.R.sh

clean:
	-rm *.pdf
	-rm /amber1/scratch/rscharpf/beaty_beadStudio/denovo.rda
	-rm /amber1/scratch/rscharpf/beaty_beadStudio/deletion.rda
	-rm /amber1/scratch/rscharpf/beaty_beadStudio/reduced_ranges.rda
