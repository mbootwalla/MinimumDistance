## A Makefile for running the cbs vignette
## must run make from a fat node
Beaty : Beaty.Rnw
	echo "Stangle(\"Beaty.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/Beaty.R/" > Beaty.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=30G,h_vmem=35G Beaty.R.sh

Crlmm : Crlmm.Rnw
	echo "Stangle(\"Crlmm.Rnw\")"| R-devel --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/Crlmm.R/" > Crlmm.R.sh
	qsub -m e -r y -cwd -l fat -l mem_free=40,h_vmem=45G Crlmm.R.sh

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
