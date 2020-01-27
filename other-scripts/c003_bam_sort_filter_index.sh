#!/bin/bash -ve
# sort, filter and index mapped read files

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=02:30:00
# sends email with maximal memory and time used by job
# current environment settings are used for the job???
#$ -V
#$ -t 1-137

taskid=${SGE_TASK_ID}

inds=`less keep_all_saxatilis.txt` # all ind IDs
ind=`echo $inds | cut -d" " -f $taskid` # get focal ind

# convert individual sam to bam, sort and index
/home/bo1awx/programs/samtools-1.2/samtools view -u -h -S $ind.sam | /home/bo1awx/programs/samtools-1.2/samtools sort - $ind
/home/bo1awx/programs/samtools-1.2/samtools index $ind.bam

# filter out secondary hits (-F 0x100 )
# filter out alignments with MAPQ smaller than 30 (-q)
# keep only targeted contigs
# then sort
/home/bo1awx/programs/samtools-1.2/samtools view -b -h -F 0x100 -q 30 $ind.bam a1234S_contig_100133 a1234S_contig_100302 a1234S_contig_100665 a1234S_contig_102211 a1234S_contig_107004 a1234S_contig_107244 a1234S_contig_108407 a1234S_contig_109408 a1234S_contig_11159 a1234S_contig_111757 a1234S_contig_112417 a1234S_contig_112696 a1234S_contig_114397 a1234S_contig_114412 a1234S_contig_11463 a1234S_contig_115731 a1234S_contig_116216 a1234S_contig_116549 a1234S_contig_117 a1234S_contig_118115 a1234S_contig_118316 a1234S_contig_118509 a1234S_contig_122784 a1234S_contig_123964 a1234S_contig_124730 a1234S_contig_125146 a1234S_contig_126474 a1234S_contig_129481 a1234S_contig_130477 a1234S_contig_13367 a1234S_contig_135492 a1234S_contig_136524 a1234S_contig_139519 a1234S_contig_140298 a1234S_contig_14442 a1234S_contig_149634 a1234S_contig_15277 a1234S_contig_153403 a1234S_contig_154806 a1234S_contig_154818 a1234S_contig_155149 a1234S_contig_15537 a1234S_contig_156218 a1234S_contig_157380 a1234S_contig_157825 a1234S_contig_1587 a1234S_contig_159909 a1234S_contig_161197 a1234S_contig_163579 a1234S_contig_166640 a1234S_contig_168338 a1234S_contig_168657 a1234S_contig_17627 a1234S_contig_177924 a1234S_contig_17881 a1234S_contig_179424 a1234S_contig_181680 a1234S_contig_186812 a1234S_contig_187672 a1234S_contig_188593 a1234S_contig_189295 a1234S_contig_189440 a1234S_contig_189481 a1234S_contig_190323 a1234S_contig_198554 a1234S_contig_198717 a1234S_contig_198904 a1234S_contig_19935 a1234S_contig_202360 a1234S_contig_203723 a1234S_contig_203788 a1234S_contig_206025 a1234S_contig_206067 a1234S_contig_206163 a1234S_contig_207417 a1234S_contig_208871 a1234S_contig_20974 a1234S_contig_215318 a1234S_contig_215451 a1234S_contig_216064 a1234S_contig_216903 a1234S_contig_217742 a1234S_contig_217743 a1234S_contig_22081 a1234S_contig_22331 a1234S_contig_224045 a1234S_contig_224569 a1234S_contig_226251 a1234S_contig_229783 a1234S_contig_230245 a1234S_contig_230253 a1234S_contig_23247 a1234S_contig_233670 a1234S_contig_23548 a1234S_contig_235780 a1234S_contig_236380 a1234S_contig_237149 a1234S_contig_239735 a1234S_contig_2401 a1234S_contig_2403 a1234S_contig_242618 a1234S_contig_242919 a1234S_contig_244198 a1234S_contig_244645 a1234S_contig_244849 a1234S_contig_245530 a1234S_contig_2460 a1234S_contig_25091 a1234S_contig_25217 a1234S_contig_255282 a1234S_contig_255534 a1234S_contig_255590 a1234S_contig_257431 a1234S_contig_257595 a1234S_contig_25804 a1234S_contig_262539 a1234S_contig_267479 a1234S_contig_269589 a1234S_contig_269725 a1234S_contig_270146 a1234S_contig_271433 a1234S_contig_271710 a1234S_contig_271817 a1234S_contig_275657 a1234S_contig_276172 a1234S_contig_279210 a1234S_contig_279742 a1234S_contig_280204 a1234S_contig_285613 a1234S_contig_289693 a1234S_contig_293233 a1234S_contig_293516 a1234S_contig_294428 a1234S_contig_296047 a1234S_contig_297497 a1234S_contig_298319 a1234S_contig_298610 a1234S_contig_300253 a1234S_contig_303270 a1234S_contig_303815 a1234S_contig_30396 a1234S_contig_305896 a1234S_contig_30614 a1234S_contig_306223 a1234S_contig_308074 a1234S_contig_309380 a1234S_contig_310095 a1234S_contig_312412 a1234S_contig_313112 a1234S_contig_3134 a1234S_contig_313564 a1234S_contig_316889 a1234S_contig_317924 a1234S_contig_318261 a1234S_contig_320 a1234S_contig_321173 a1234S_contig_323421 a1234S_contig_324661 a1234S_contig_328131 a1234S_contig_3290 a1234S_contig_331154 a1234S_contig_33142 a1234S_contig_333081 a1234S_contig_333949 a1234S_contig_33460 a1234S_contig_340045 a1234S_contig_34235 a1234S_contig_343664 a1234S_contig_34375 a1234S_contig_3438 a1234S_contig_34627 a1234S_contig_348565 a1234S_contig_348615 a1234S_contig_35062 a1234S_contig_350735 a1234S_contig_353340 a1234S_contig_355000 a1234S_contig_357034 a1234S_contig_35946 a1234S_contig_361678 a1234S_contig_363598 a1234S_contig_363698 a1234S_contig_364953 a1234S_contig_37210 a1234S_contig_372596 a1234S_contig_372640 a1234S_contig_3739 a1234S_contig_3790 a1234S_contig_379419 a1234S_contig_379710 a1234S_contig_37997 a1234S_contig_380025 a1234S_contig_381939 a1234S_contig_381993 a1234S_contig_385767 a1234S_contig_389608 a1234S_contig_391508 a1234S_contig_3965 a1234S_contig_397649 a1234S_contig_397663 a1234S_contig_398202 a1234S_contig_39851 a1234S_contig_40080 a1234S_contig_402912 a1234S_contig_403851 a1234S_contig_408777 a1234S_contig_40975 a1234S_contig_411529 a1234S_contig_414919 a1234S_contig_415943 a1234S_contig_416537 a1234S_contig_421256 a1234S_contig_421633 a1234S_contig_428361 a1234S_contig_42855 a1234S_contig_434479 a1234S_contig_437469 a1234S_contig_437506 a1234S_contig_440349 a1234S_contig_441082 a1234S_contig_44162 a1234S_contig_441934 a1234S_contig_446245 a1234S_contig_44645 a1234S_contig_451636 a1234S_contig_455668 a1234S_contig_457402 a1234S_contig_46447 a1234S_contig_464797 a1234S_contig_465846 a1234S_contig_469023 a1234S_contig_469414 a1234S_contig_474549 a1234S_contig_476390 a1234S_contig_478500 a1234S_contig_478844 a1234S_contig_48392 a1234S_contig_484050 a1234S_contig_486679 a1234S_contig_486721 a1234S_contig_486797 a1234S_contig_4893 a1234S_contig_492559 a1234S_contig_492877 a1234S_contig_4980 a1234S_contig_49954 a1234S_contig_50160 a1234S_contig_504691 a1234S_contig_505817 a1234S_contig_506140 a1234S_contig_509081 a1234S_contig_5138 a1234S_contig_51969 a1234S_contig_521140 a1234S_contig_524279 a1234S_contig_52456 a1234S_contig_524995 a1234S_contig_52534 a1234S_contig_525714 a1234S_contig_525942 a1234S_contig_530709 a1234S_contig_532610 a1234S_contig_535448 a1234S_contig_544335 a1234S_contig_54577 a1234S_contig_550359 a1234S_contig_55097 a1234S_contig_554703 a1234S_contig_554831 a1234S_contig_557256 a1234S_contig_560678 a1234S_contig_562681 a1234S_contig_56375 a1234S_contig_56605 a1234S_contig_56665 a1234S_contig_569763 a1234S_contig_57137 a1234S_contig_573958 a1234S_contig_582801 a1234S_contig_587601 a1234S_contig_588329 a1234S_contig_590061 a1234S_contig_590933 a1234S_contig_592488 a1234S_contig_59363 a1234S_contig_594424 a1234S_contig_59604 a1234S_contig_596699 a1234S_contig_598398 a1234S_contig_601374 a1234S_contig_614450 a1234S_contig_61475 a1234S_contig_616193 a1234S_contig_624013 a1234S_contig_62621 a1234S_contig_628668 a1234S_contig_6372 a1234S_contig_637899 a1234S_contig_64055 a1234S_contig_648296 a1234S_contig_648767 a1234S_contig_653267 a1234S_contig_661592 a1234S_contig_662892 a1234S_contig_664702 a1234S_contig_665647 a1234S_contig_6660 a1234S_contig_666366 a1234S_contig_666513 a1234S_contig_670067 a1234S_contig_67106 a1234S_contig_674122 a1234S_contig_67604 a1234S_contig_676623 a1234S_contig_68357 a1234S_contig_683715 a1234S_contig_684481 a1234S_contig_689043 a1234S_contig_691779 a1234S_contig_70718 a1234S_contig_70869 a1234S_contig_716489 a1234S_contig_716637 a1234S_contig_72023 a1234S_contig_7208 a1234S_contig_72221 a1234S_contig_74805 a1234S_contig_7567 a1234S_contig_76500 a1234S_contig_79560 a1234S_contig_80053 a1234S_contig_847 a1234S_contig_85770 a1234S_contig_86094 a1234S_contig_87117 a1234S_contig_87190 a1234S_contig_88656 a1234S_contig_92401 a1234S_contig_93974 a1234S_contig_94541 a1234S_contig_95561 a1234S_contig_95778 a1234S_contig_96046 a1234S_contig_98628 | /home/bo1awx/programs/samtools-1.2/samtools sort - $ind-filtered

# index filtered & sorted file, and remove unfiltered versions
/home/bo1awx/programs/samtools-1.2/samtools index $ind-filtered.bam
rm $ind.bam $ind.bam.bai