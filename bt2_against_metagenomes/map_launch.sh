
curDir=$(realpath $(dirname $0));



odir="/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/prevalence2/";
prefix="/shares/CIBIO-Storage/CM/scratch/data/meta"
base=$prefix;
##############################################
 
mode=$1;

if [ $mode == "U" ]; then
	extension='.fastq'
	extraction_cmd='cat'
fi;

if [ $mode == "bzip" ]; then
	extension='.fastq.bz2'
	extraction_cmd='bzcat'
fi;
	
if [ $mode == "gzip" ]; then
	extension='.fastq.gz'
	extraction_cmd='zcat'
fi;

if [ $mode == "fqgz" ]; then
	extension='.fq.gz'
	extraction_cmd='zcat'
fi;

etp='HMP_2012
LiuB_2012
BengtssonPalmeJ_2015
PiperHG_2016
ContevilleLC_2019
GopalakrishnanV_2018
RosaBA_2018
KieserS_2018
LassalleF_2017
SankaranarayananK_2015
RampelliS_2015
MatsonV_2018
SmitsSA_2017
WingleeK_2017
LomanNJ_2013
CM_ghana
OlmMR_2017
LawrenceA_2015
CM_ethiopia
GeversD_2014
Heitz-BuschartA_2016
LiSS_2016
LokmerA_2019
Obregon-TitoAJ_2015
CM_rescignocrc
GuptaA_2019
WampachL_2018
VoigtAY_2015
CM_tanzania
RaymondF_2016
YeZ_2018
KorpelaK_2016
ThomasAM_2019_c
HanniganGD_2017
LoombaR_2017
LouisS_2016
CM_lilt
ChengpingW_2017
IjazUZ_2017
ChuDM_2017
CM_madagascar
DhakanDB_2019
LiuW_2016
VogtmannE_2016
CM_guinea
HeQ_2017
KosticAD_2015
WirbelJ_2018
YuJ_2015
KarlssonFH_2013
FengQ_2015
ZellerG_2014
WenC_2017
LiJ_2017
HansenLBS_2018
VincentC_2016
QinN_2014
CM_caritro
XieH_2016
LiJ_2014
HallAB_2017
CosteaPI_2017
PehrssonE_2016
YassourM_2018
LeChatelierE_2013
QinJ_2012
JieZ_2017
NielsenHB_2014
BackhedF_2015
SchirmerM_2016
CM_caritro_twins
CM_tanzania2
CM_sardegna
CM_guinea2
CM_ghana2'

for et in $etp; do
	for utp in $(find ${prefix}/${et}/reads -maxdepth 1 -mindepth 1 -type d); do
		e+=($utp);
	done;
done;

while true; do
	
	bt=0;
	b_short=$(qstat -u moreno.zolfo | grep short | grep "VDBM" | wc -l);
	b_common=$(qstat -u moreno.zolfo | grep common | grep "VDBM" | wc -l);
	b_cibio=$(qstat -u moreno.zolfo | grep CIBIO_ | grep "VDBM" | wc -l);
	b_cibiocm=$(qstat -u moreno.zolfo | grep CIBIOCM | grep "VDBM" | wc -l);

	
	for folder in ${e[@]}; do
		
		sample=$(basename $folder);
		dataset=$(basename $(dirname $(dirname $folder)));
		mkdir -p ${odir}/${dataset}/
		odirE=${odir}/${dataset}/

		if [ ! -f ${odir}/${dataset}/vdbm__${dataset}__${sample}.bam  ] && [ ! -f ${odir}/${dataset}/vdbm__${dataset}__${sample}.wk ]; then 		
			if [ $b_short -lt 30 ]; then 
				#echo "L" ${odir}/${dataset}/vdbm__${dataset}__${sample}.wk
				qsub -q short_cpuQ -v prefix=\"${base}\",outDir=\"${odirE}\",dataset=\"${dataset}\",sample=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"16\" -N VDBM_${dataset}_${sample} -l select=1:ncpus=16 ${curDir}/map_bt2.sh;
				touch ${odir}/${dataset}/vdbm__${dataset}__${sample}.wk
				b_short=$((b_short+1));
				bt=$((bt+1));
			elif [ $b_common -lt 50 ]; then 
				#echo "L" ${odir}/${dataset}/vdbm__${dataset}__${sample}.wk
				qsub -q common_cpuQ -v prefix=\"${base}\",outDir=\"${odirE}\",dataset=\"${dataset}\",sample=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"8\" -N VDBM_${dataset}_${sample} -l select=1:ncpus=8 ${curDir}/map_bt2.sh;
				touch ${odir}/${dataset}/vdbm__${dataset}__${sample}.wk
				b_common=$((b_common+1));
				bt=$((bt+1)); 
			elif [ $b_cibio -lt 50 ]; then 
				#echo "L" ${odir}/${dataset}/vdbm__${dataset}__${sample}.wk
				qsub -q CIBIO_cpuQ -v prefix=\"${base}\",outDir=\"${odirE}\",dataset=\"${dataset}\",sample=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"8\" -N VDBM_${dataset}_${sample} -l select=1:ncpus=8 ${curDir}/map_bt2.sh;
				touch ${odir}/${dataset}/vdbm__${dataset}__${sample}.wk
				b_cibio=$((b_cibio+1));
				bt=$((bt+1));
			elif [ $b_cibiocm -lt 50 ]; then 
				#echo "L" ${odir}/${dataset}/vdbm__${dataset}__${sample}.wk
				qsub -q CIBIOCM_cpuQ -v prefix=\"${base}\",outDir=\"${odirE}\",dataset=\"${dataset}\",sample=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"8\" -N VDBM_${dataset}_${sample} -l select=1:ncpus=8 ${curDir}/map_bt2.sh;
				touch ${odir}/${dataset}/vdbm__${dataset}__${sample}.wk
				b_cibiocm=$((b_cibiocm+1));
				bt=$((bt+1));
			fi;
		fi;
	done;

	a=$(qstat -u moreno.zolfo | grep "VDBM" | grep " R " | wc -l);
	q=$(qstat -u moreno.zolfo | grep "VDBM" | grep " Q " | wc -l);
	ete=${#e[@]}
	ft="ND" #$(find /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/vqc/ -name "*.vqc.txt" | wc -l);
	
	tg="ND" #$(($ete-$ft));
	echo $(date) "| " $bt "jobs launched, " $a "running ", $q "queued. See you in one minute." ;
	sleep 60
	b=0;
done;
