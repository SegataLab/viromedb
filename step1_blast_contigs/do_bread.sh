#!/bin/bash

DATA="/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/contigs_tg/original/";

export curDir=$(realpath $(dirname $0));
export DATA;

vtp="c7000-1__100__LiangG_2020__D3334__D3334
c7001-1__100__LiangG_2020__D100__D100
c7002-1__100__LiangG_2020__D101__D101
c7003-1__100__LiangG_2020__D105__D105
c7004-1__64__LiangG_2020__D126__D126
c7005-1__100__LiangG_2020__D133__D133
c7006-1__57__LiangG_2020__D145__D145
c7007-1__100__LiangG_2020__D196__D196
c7008-1__87__LiangG_2020__D215__D215
c7009-1__100__LiangG_2020__D230__D230
c7010-1__89__LiangG_2020__D233__D233
c7011-1__100__LiangG_2020__D234__D234
c7012-1__100__LiangG_2020__D245__D245
c7013-1__100__LiangG_2020__D261__D261
c7014-1__75__LiangG_2020__D275__D275
c7015-1__100__LiangG_2020__D280__D280
c7016-1__98__LiangG_2020__D283__D283
c7017-1__100__LiangG_2020__D287__D287
c7018-1__73__LiangG_2020__D306__D306
c7019-1__100__LiangG_2020__D3160__D3160
c7020-1__73__LiangG_2020__D3162__D3162
c7021-1__83__LiangG_2020__D3164__D3164
c7022-1__100__LiangG_2020__D3170__D3170
c7023-1__100__LiangG_2020__D3175__D3175
c7024-1__100__LiangG_2020__D3199__D3199
c7025-1__100__LiangG_2020__D3251__D3251
c7026-1__100__LiangG_2020__D3255__D3255
c7027-1__51__LiangG_2020__D3262__D3262
c7028-1__100__LiangG_2020__D3264__D3264
c7029-1__100__LiangG_2020__D3265__D3265
c7030-1__57__LiangG_2020__D3267__D3267
c7031-1__66__LiangG_2020__D3273__D3273
c7032-1__100__LiangG_2020__D3280__D3280
c7033-1__100__LiangG_2020__D3281__D3281
c7034-1__60__LiangG_2020__D3451__D3451
c7035-1__92__LiangG_2020__D3454__D3454
c7036-1__100__LiangG_2020__D348__D348
c7038-1__100__LiangG_2020__D358__D358
c7039-1__100__LiangG_2020__D3594__D3594
c7040-1__100__LiangG_2020__D3681__D3681
c7042-1__100__LiangG_2020__D370__D370
c7044-1__76__LiangG_2020__D3764__D3764
c7045-1__100__LiangG_2020__D3854__D3854
c7046-1__100__LiangG_2020__D386__D386
c7047-1__100__LiangG_2020__D392__D392
c7048-1__100__LiangG_2020__D3941__D3941
c7049-1__72__LiangG_2020__D3951__D3951
c7050-1__84__LiangG_2020__D3964__D3964
c7051-1__100__LiangG_2020__D4001__D4001
c7052-1__100__LiangG_2020__D4014__D4014
c7053-1__51__LiangG_2020__D4284__D4284
c7054-1__100__LiangG_2020__DB5__DB5
c7055-1__100__LiangG_2020__DD2__DD2
c7056-1__58__LiangG_2020__DDH__DDH
c7057-1__100__LiangG_2020__DDN2__DDN2
c7058-1__100__LiangG_2020__yD1015__yD1015
c7059-1__100__LiangG_2020__yD1171__yD1171
c7060-1__100__LiangG_2020__yD1306__yD1306
c7061-1__51__LiangG_2020__yD1440__yD1440
c7062-1__100__LiangG_2020__yD1582__yD1582
c7063-1__100__LiangG_2020__yD171__yD171
c7064-1__100__LiangG_2020__yD1936__yD1936
c7065-1__100__LiangG_2020__yD378__yD378
c7066-1__65__LiangG_2020__yD393__yD393
c7067-1__100__LiangG_2020__yD408__yD408
c7068-1__100__LiangG_2020__yD519__yD519
c7069-1__100__LiangG_2020__yD526__yD526
c7073-1__64__LiangG_2020__D4281__D4281
c7099-1__60__LiangG_2020__D3944__D3944";

vtp="c7037-1__60__LiangG_2020__D3521__D3521
c7041-1__67__LiangG_2020__D3684__D3684
c7043-1__100__LiangG_2020__D3734__D3734"

vtp="c8000-1__100__RefSeq__RefSeq__RefSeq"

parallel -j 22 --env curDir --env DATA 'i={}; ${curDir}/bread.py ${DATA}/${i}.ncbi.blast#NCBI80k#80#1000 ${DATA}/${i}.vir91.blast#RefSeq#80#500 ${DATA}/${i}.sgbs.blast#SGBS#80#1000 > ./${i}.csv' ::: $vtp;


#c7037-1__60__LiangG_2020__D3521__D3521
#c7041-1__67__LiangG_2020__D3684__D3684
#c7043-1__100__LiangG_2020__D3734__D3734
