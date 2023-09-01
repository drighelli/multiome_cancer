ALL: HLTY CNCR

wrkrs = 8

HLTY: REP1_HLTY REP2_HLTY REP3_HLTY

CNCR: REP1_CNCR REP2_CNCR 

REP1_CNCR:
	cellranger-arc count --id=$@ \
	--reference=/mnt/ganimede/data/references/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
	--libraries=REP1_CNCR.csv \
	--localcores=$(wrkrs) \
	--localmem=40
	
REP2_CNCR:
	cellranger-arc count --id=$@ \
	--reference=/mnt/ganimede/data/references/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
	--libraries=REP2_CNCR.csv \
	--localcores=$(wrkrs) \
	--localmem=40

REP3_CNCR:
	cellranger-arc count --id=$@ \
	--reference=/mnt/ganimede/data/references/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
	--libraries=REP3_CNCR.csv \
	--localcores=$(wrkrs) \
	--localmem=40
	
REP1_HLTY:
	cellranger-arc count --id=$@ \
	--reference=/mnt/ganimede/data/references/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
	--libraries=REP1_HLTY.csv \
	--localcores=$(wrkrs) \
	--localmem=40
	
REP2_HLTY:
	cellranger-arc count --id=$@ \
	--reference=/mnt/ganimede/data/references/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
	--libraries=REP2_HLTY.csv \
	--localcores=$(wrkrs) \
	--localmem=40

REP3_HLTY:
	cellranger-arc count --id=$@ \
	--reference=/mnt/ganimede/data/references/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
	--libraries=REP3_HLTY.csv \
	--localcores=$(wrkrs) \
	--localmem=40
