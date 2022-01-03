#! /bin/sh

# sorhum reference
#for SAMPLEID in 22260_51 22261_71 22264_51 22264_71 22264_91 22269_51 22274_51 22274_71 #22275_71 22275_91 22276_51 22277_51 22277_71 22277_91 22278_11 22279_91 22282_71 22283_51 #22283_71 22283_91 22295_11 22297_91 22298_91 22300_71 22302_71 22302_91 22304_71 22305_91 #22306_71 22318_71 22320_71 22320_91 22321_91 22322_51 22322_71 22322_91 22323_51 #22324_91_ 22324_11 22325_71 22325_91 22330_71 22331_71 22332_71 22332_91 22332_11 #22334_91 22335_11 22336_71 22336_91 22337_71 22338_91 22340_11 22341_51 22341_71 22342_91 #22343_91 22345_91 22346_91 22347_91 22348_71 22349_91 22352_71 22352_91 22355_71 22355_91

# spinach reference
for SAMPLEID in 22260_51 22261_71 22262_91 22264_51 22264_71 22264_91 22266_91 22267_111 22269_51 22274_51 22274_71 22275_71 22275_91 22277_51 22277_71 22277_91 22278_111 22279_91 22282_71 22283_51 22283_71 22283_91 22295_111 22297_91 22298_91 22300_71 22302_71 22302_91 22304_71 22305_71 22305_91 22306_71 22318_71 22320_71 22320_91 22321_91 22322_51 22322_71 22322_91 22323_51 22324_91 22324_111 22325_71 22325_91 22328_71 22330_71 22331_71 22332_71 22332_91 22332_111 22334_91 22335_111 22336_71 22336_91 22337_71 22338_91 22340_111 22341_51 22341_71 22342_91 22343_91 22345_91 22346_91 22347_91 22348_71 22349_91 22352_71 22352_91 22355_71 22355_91

do

~/groupdata/jre1.8.0_51/bin/java -jar Chloe.jar $SAMPLEID/contigs.fa $SAMPLEID/So.coords $SAMPLEID/$SAMPLEID\_So.chloe

done