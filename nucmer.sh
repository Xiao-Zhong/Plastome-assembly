#! /bin/sh

for SAMPLEID in 22260 22261 22262 22263 22264 22265 22266 22267 22268 22269 22270 22271 22272 22273 22274 22275 22276 22277 22278 22279 22280 22281 22282 22283 22284 22285 22286 22287 22288 22289 22290 22291 22292 22293 22294 22295 22296 22297 22298 22299 22300 22301 22302 22303 22304 22305 22306 22307 22308 22309 22310 22311 22312 22313 22314 22315 22316 22317 22318 22319 22320 22321 22322 22323 22324 22325 22326 22327 22328 22329 22330 22331 22332 22333 22334 22335 22336 22337 22338 22339 22340 22341 22342 22343 22344 22345 22346 22347 22348 22349 22350 22351 22352 22353 22354 22355

do

/usr/local/MUMmer3.23/nucmer --prefix=$SAMPLEID\_51/Sb EF115542.fasta $SAMPLEID\_51/contigs.fa
/usr/local/MUMmer3.23/show-coords -rcl $SAMPLEID\_51/Sb.delta > $SAMPLEID\_51/Sb.coords
/usr/local/MUMmer3.23/nucmer --prefix=$SAMPLEID\_71/Sb EF115542.fasta $SAMPLEID\_71/contigs.fa
/usr/local/MUMmer3.23/show-coords -rcl $SAMPLEID\_71/Sb.delta > $SAMPLEID\_71/Sb.coords
/usr/local/MUMmer3.23/nucmer --prefix=$SAMPLEID\_91/Sb EF115542.fasta $SAMPLEID\_91/contigs.fa
/usr/local/MUMmer3.23/show-coords -rcl $SAMPLEID\_91/Sb.delta > $SAMPLEID\_91/Sb.coords
/usr/local/MUMmer3.23/nucmer --prefix=$SAMPLEID\_111/Sb EF115542.fasta $SAMPLEID\_111/contigs.fa
/usr/local/MUMmer3.23/show-coords -rcl $SAMPLEID\_111/Sb.delta > $SAMPLEID\_111/Sb.coords

/usr/local/MUMmer3.23/nucmer --prefix=$SAMPLEID\_51/So NC_002202.fasta $SAMPLEID\_51/contigs.fa
/usr/local/MUMmer3.23/show-coords -rcl $SAMPLEID\_51/So.delta > $SAMPLEID\_51/So.coords
/usr/local/MUMmer3.23/nucmer --prefix=$SAMPLEID\_71/So NC_002202.fasta $SAMPLEID\_71/contigs.fa
/usr/local/MUMmer3.23/show-coords -rcl $SAMPLEID\_71/So.delta > $SAMPLEID\_71/So.coords
/usr/local/MUMmer3.23/nucmer --prefix=$SAMPLEID\_91/So NC_002202.fasta $SAMPLEID\_91/contigs.fa
/usr/local/MUMmer3.23/show-coords -rcl $SAMPLEID\_91/So.delta > $SAMPLEID\_91/So.coords
/usr/local/MUMmer3.23/nucmer --prefix=$SAMPLEID\_111/So NC_002202.fasta $SAMPLEID\_111/contigs.fa
/usr/local/MUMmer3.23/show-coords -rcl $SAMPLEID\_111/So.delta > $SAMPLEID\_111/So.coords

done