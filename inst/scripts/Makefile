## A Makefile for running the cbs vignette
## must run make from a fat node
Beaty : Beaty.Rnw
	echo "Stangle(\"Beaty.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/Beaty.R/" > Beaty.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=30G,h_vmem=35G Beaty.R.sh

segmentation: segmentation.Rnw
	echo "Stangle(\"segmentation.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/segmentation.R/" > segmentation.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=5G,h_vmem=12G segmentation.R.sh

readBS: ReadBeadStudio.Rnw	
	echo "Stangle(\"ReadBeadStudio.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/ReadBeadStudio.R/" > ReadBeadStudio.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=5G,h_vmem=12G ReadBeadStudio.R.sh

MinDistance: MinDistance.Rnw
	echo "Stangle(\"MinDistance.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/MinDistance.R/" > MinDistance.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=5G,h_vmem=12G MinDistance.R.sh

PennCNV: PennCNV.Rnw
	echo "Stangle(\"PennCNV.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/PennCNV.R/" > PennCNV.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=5G,h_vmem=12G PennCNV.R.sh

figs: figs.Rnw
	echo "Stangle(\"figs.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/figs.R/" > figs.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=1G,h_vmem=2G figs.R.sh

DenovoRanges: DenovoRangesBySize.Rnw
	echo "Stangle(\"DenovoRangesBySize.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/DenovoRangesBySize.R/" > DenovoRangesBySize.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=25G,h_vmem=30G DenovoRangesBySize.R.sh


clean:
	-rm *.pdf
	-rm /amber1/scratch/rscharpf/beaty_beadStudio/denovo.rda
	-rm /amber1/scratch/rscharpf/beaty_beadStudio/deletion.rda
	-rm /amber1/scratch/rscharpf/beaty_beadStudio/reduced_ranges.rda
