#!/bin/bash

#############################################################################################################
#                                                                                                           #
# JolyTree: fast distance-based phylogenetic inference from unaligned genome sequences                      #
#                                                                                                           #
  COPYRIGHT="Copyright (C) 2017-2021 Institut Pasteur"                                                      #
#                                                                                                           #
# This program  is free software:  you can  redistribute it  and/or modify it  under the terms  of the GNU  #
# General Public License as published by the Free Software Foundation, either version 3 of the License, or  #
# (at your option) any later version.                                                                       #
#                                                                                                           #
# This program is distributed in the hope that it will be useful,  but WITHOUT ANY WARRANTY;  without even  #
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public  #
# License for more details.                                                                                 #
#                                                                                                           #
# You should have received a copy of the  GNU General Public License along with this program.  If not, see  #
# <http://www.gnu.org/licenses/>.                                                                           #
#                                                                                                           #
# Contact:                                                                                                  #
#  Alexis Criscuolo                                                            alexis.criscuolo@pasteur.fr  #
#  Genome Informatics & Phylogenetics (GIPhy)                                             giphy.pasteur.fr  #
#  Bioinformatics and Biostatistics Hub                                 research.pasteur.fr/team/hub-giphy  #
#  USR 3756 IP CNRS                          research.pasteur.fr/team/bioinformatics-and-biostatistics-hub  #
#  Dpt. Biologie Computationnelle                     research.pasteur.fr/department/computational-biology  #
#  Institut Pasteur, Paris, FRANCE                                                     research.pasteur.fr  #
#                                                                                                           #
#############################################################################################################

#############################################################################################################
#                                                                                                           #
# ============                                                                                              #
# = VERSIONS =                                                                                              #
# ============                                                                                              #
#                                                                                                           #
  VERSION=2.1.211019ac                                                                                      #
# + commenting line 576, as some linux distribution incorrectly interpret  "trap [arg] signal_spec"  with   #
#   empty arg                                                                                               #
# + adding some conditions to deal with some FastME crashes observed when inferring large trees             #
# + new option -x to prevent 
#                                                                                                           #
# VERSION=2.0.190926ac                                                                                      #
# + new F81/EI transformation formula using gamma shape parameter (option -a = 1.5 by default)              #
# + option -f to to use the 4 nucleotide frequencies in F81/EI transformation; by default, to deal with     #
#   multiple contig files, JolyTree sets f(A)=f(T)=0.5*(A+T)/(A+C+G+T) and f(C)=f(G)=0.5*(C+G)/(A+C+G+T)    #
# + option -s modified                                                                                      #
#                                                                                                           #
# VERSION=1.2.190729ac                                                                                      #
# + 3x faster minhash estimations on multiple threads                                                       #
# + ratchet-based BME tree search based on multiple threads                                                 #
# + temporary files deleted when JolyTree is interrupted                                                    #
#                                                                                                           #
# VERSION=1.1.181205ac                                                                                      #
# + option -q to set desired probability of observing a random k-mer                                        #
# + ability to be run on clusters managed by SLURM                                                          #
#                                                                                                           #
# VERSION=1.0.180115ac                                                                                      #
# + option -n to only estimate evolutionary distances                                                       #
# + option -r to set the number of iterations when performing the ratchet-based BME tree search             #
# + important bug fixed when sorting the input file names                                                   #
# + reimplementation of the F81 distance estimation                                                         #
#                                                                                                           #
# VERSION=0.8.171207ac                                                                                      #
# + k-mer size could be set by the user                                                                     #
# + bug fixed for manual tbl estimation                                                                     #
# + .fas allowed inside the input directory                                                                 #
#                                                                                                           #
# VERSION=0.7.170919ac                                                                                      #
# + tree output file suffix is now .nwk                                                                     #
#                                                                                                           #
# VERSION=0.6.170728ac                                                                                      #
# + implements the F81  transformation suggested  by Tamura & Kumar (2002)  in order to deal with putative  #
#   heterogeneous substitution pattern among lineages                                                       #
#                                                                                                           #
# VERSION=0.5.170727ac                                                                                      #
# + automatic estimation of the k-mer size                                                                  #
#                                                                                                           #
# VERSION=0.4.170726ac                                                                                      #
# + precomputed pairwise p-distances could be used (option -d)                                              #
# + no limit with the length of the input FASTA filenames                                                   #
#                                                                                                           #
# VERSION=0.3.170724ac                                                                                      #
# + implements the F81 transformation when at least one p-distance is larger than a specified cutoff        #
#                                                                                                           #
# VERSION=0.2.170721ac                                                                                      #
# + uses FastME 2.1.5.1                                                                                     #
#                                                                                                           #
#############################################################################################################

#############################################################################################################
#                                                                                                           #
# ============                                                                                              #
# = DOC      =                                                                                              #
# ============                                                                                              #
#                                                                                                           #
  if [ "$1" = "-?" ] || [ "$1" = "-h" ] || [ $# -le 1 ]                                                     #
  then                                                                                                      #
    echo -e "\n\033[1m JolyTree v$VERSION                         $COPYRIGHT\033[0m";                       #
    cat << EOF                                                       

 Criscuolo A  (2019)  A fast alignment-free  bioinformatics procedure to  infer accurate 
 distance-based phylogenetic trees from genome assemblies. RIO. doi:10.3897/rio.5.e36178 

 Criscuolo A  (2020)  On the transformation of  MinHash-based uncorrected distances into 
 proper   evolutionary    distances   for    phylogenetic   inference.    F1000Research.
 doi:10.12688/f1000research.26930.1

 USAGE:  JolyTree.sh  -i <directory>  -b <basename>  [options]

 OPTIONS:

    -i <directory>  directory name containing  FASTA-formatted contig files;  only files
                    ending with .fa, .fna, .fas or .fasta will be considered (mandatory)
    -b <basename>   basename of every written output file (mandatory)
    -q <real>       probability of observing a random k-mer (default: 0.000000001)
    -k <int>        k-mer size (default: estimated from the largest genome size with the
                    probability set by option -q)
    -s <real|int>   sketch size estimated as the  proportion (up to 1.00) of the average 
                    genome size;  the sketch size (integer > 2) can also be directly set 
                    using this option (default: 0.25)
    -c <real>       if at least one of the estimated p-distances is above this specified
                    cutoff, then a F81/EI correction is performed (default: 0.1)
    -a <real>       F81/EI transformation gamma shape parameter (default: 1.5)
    -f              using  the  four nucleotide  frequencies in  F81/EI  transformations  
                    (default:  to  deal  with  multiple  contig  files,   only  the  two 
                    frequencies A+T and C+G are used)
    -n              no BME tree inference (only pairwise distance estimates)
    -r <int>        number of steps  when performing the  ratchet-based  BME tree search
                    (default: 100)
    -x              no branch support
    -t <int>        number of threads (default: 2)

EOF
    if [ "$1" = "-?" ] || [ "$1" = "-h" ] || [ $# -le 0 ]; then exit 0 ; fi                                 #
    exit 1 ;                                                                                                #
  fi                                                                                                        #
#                                                                                                           #
#############################################################################################################

  
#############################################################################################################
#                                                                                                           #
# ================                                                                                          #
# = INSTALLATION =                                                                                          #
# ================                                                                                          #
#                                                                                                           #
# [1] REQUIREMENTS =======================================================================================  #
#  JolyTree depends on Mash, gawk, FastME and REQ (see below),  each with a minimum version required. You   #
#  should have them installed on your computer prior to using JolyTree.  Make sure that each is installed   #
#  on your $PATH variable, or specify below the full path to each of them.                                  #
#                                                                                                           #
# -- Mash: fast pairwise p-distance estimation --------------------------------------------------------     #
#    VERSION >= 1.0.2                                                                                       #
#    src: github.com/marbl/Mash                                                                             #
#    Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM (2016) Mash: fast     #
#      genome  and  metagenome  distance  estimation  using  MinHash.   Genome  Biology,  17:132.  doi:     #
#      10.1186/s13059-016-0997-x                                                                            #
#                                                                 ################################################
                                                                  ################################################
  MASH=mash;                                                      ## <=== WRITE HERE THE PATH TO THE MASH       ##
                                                                  ##      BINARY (VERSION 1.0.2 MINIMUM)        ##
                                                                  ################################################
                                                                  ################################################
#                                                                                                           #
# -- gawk: fast text file processing ------------------------------------------------------------------     #
#    VERSION >= 4.1.0                                                                                       #
#    src: ftp.gnu.org/gnu/gawk                                                                              #
#    Robbins AD  (2018)  GAWK:  Effective AWK Programming -- A User’s Guide  for GNU Awk  (Edition 4.2)     #
#      www.gnu.org/software/gawk/manual                                                                     #
#                                                                 ################################################
                                                                  ################################################
  GAWK=gawk;                                                      ## <=== WRITE HERE THE PATH TO THE GAWK       ##
                                                                  ##      BINARY (VERSION 4.1.0 MINIMUM)        ##
                                                                  ################################################
                                                                  ################################################
#                                                                                                           #
# -- FastME: fast distance-based phylogenetic tree inference ------------------------------------------     #
#    VERSION >= 2.1.5.1                                                                                     #
#    src: gite.lirmm.fr/atgc/FastME/                                                                        #
#    Lefort V, Desper R, Gascuel O  (2015)  FastME 2.0:  a comprehensive, accurate,  and fast distance-     #
#      based  phylogeny  inference  program.  Molecular Biology and Evolution,  32(10):2798–2800.  doi:     #
#      10.1093/molbev/msv150                                                                                #
#                                                                 ################################################
                                                                  ################################################
  FASTME=fastme;                                                  ## <=== WRITE HERE THE PATH TO THE FASTME     ##
                                                                  ##      BINARY (VERSION 2.1.5.1 MINIMUM)      ##
                                                                  ################################################
                                                                  ################################################
#                                                                                                           #
# -- REQ: fast computation of the rates of elementary quartets ----------------------------------------     #
#    VERSION >= 1.2                                                                                         #
#    src: gitlab.pasteur.fr/GIPhy/REQ                                                                       #
#    Guenoche A, Garreta H (2001) Can we have confidence in a tree representation. In: Gascuel O, Sagot     #
#      MF (eds) Computational Biology.  Lecture Notes in Computer Science, vol 2066.  Springer, Berlin,     #
#      Heidelberg. doi:10.1007/3-540-45727-5_5                                                              #
#                                                                 ################################################
                                                                  ################################################
  REQ=REQ;                                                        ## <=== WRITE HERE THE PATH TO THE REQ        ##
                                                                  ##      BINARY (VERSION 1.2 MINIMUM)          ##
                                                                  ################################################
                                                                  ################################################
#                                                                                                           #
#                                                                                                           #
#                                                                                                           #
#                                                                                                           #
# [2] EXECUTE PERMISSION =================================================================================  #
#  In order to run JolyTree, give the execute permission on the script JolyTree.sh:                         #
#    chmod +x JolyTree.sh                                                                                   #
#                                                                                                           #
#                                                                                                           #
#                                                                                                           #
# [3] NOTES ON THE USE OF JOLYTREE WITH SLURM (slurm.schedmd.com) ========================================  #
#  By default, JolyTree is able to perform the pairwise p-distance estimate step on multiple threads (the   #
#  option -t  allows the  number of  threads  to  be specified).  The corresponding  pieces of  codes are   #
#  therefore executed concurrently via the following standard 'launcher':                                   #
#                                                                                                           #
   EXEC="sh -c";                                                                                            #
#                                                                                                           #
#  It is therefore possible to use JolyTree on a computer that allows multiple threads to be executed. It   #
#  is  also possible  to launch  JolyTree on  multiple  threads  on a  cluster managed  by Slurm  via the   #
#  following command line models (with t = number of threads):                                              #
#    srun   <Slurm options> -c $t  ./JolyTree.sh  <JolyTree options>  -t $t                                 #
#    sbatch <Slurm options> -c $t  ./JolyTree.sh  <JolyTree options>  -t $t                                 #
#  Moreover, it is also possible to launch JolyTree on  multiple cores on a cluster managed by Slurm. For   #
#  this particular case, you should first uncomment the following line:                                     #
#                                                                                                           #
#  EXEC="srun -n 1 -N 1 -Q $EXEC";                                                                          #
#                                                                                                           #
#  and launch JolyTree via the following command line models (with t = number of cores):                    #
#    salloc <Slurm options> -n $t  ./JolyTree.sh  <JolyTree options>  -t $t                                 #
#    sbatch <Slurm options> -n $t  ./JolyTree.sh  <JolyTree options>  -t $t                                 #
#                                                                                                           #
#############################################################################################################

  
#############################################################################################################
#############################################################################################################
#### INITIALIZING PARAMETERS AND READING OPTIONS                                                         ####
#############################################################################################################
#############################################################################################################

if [ ! $(command -v $MASH) ];   then echo "$MASH not found"   >&2 ; exit 1 ; fi
if [ ! $(command -v $GAWK) ];   then echo "$GAWK not found"   >&2 ; exit 1 ; fi
if [ ! $(command -v $FASTME) ]; then echo "$FASTME not found" >&2 ; exit 1 ; fi
if [ ! $(command -v $REQ) ];    then echo "$REQ not found"    >&2 ; exit 1 ; fi

export LC_ALL=C ;

DATADIR="N.O_D.I.R";            # -i (mandatory)
BASEFILE="N.O_B.A.S.E.F.I.L.E"; # -b (mandatory)

SKETCH=0.25;                    # -s (0.25)
Q=0.000000001;                  # -q (0.000000001)
K=0;                            # -k (auto from -q)
CUTOFF=0.1;                     # -c (0.1)
ALPHA=1.5;                      # -a (1.5)
NFQ=2;                          # -f (none)

INFERTREE=true;                 # -n (none)
RATCHET=100;                    # -r (100)
RATCHET_LIMIT=200;              #    (static)

BRANCH_SUPPORT=true;            # -x (none)

NPROC=2;                        # -t (2)
CHUNK=20;                       # -h (20)
WAITIME=0.5;                    #    (auto from -t)

while getopts :i:b:s:q:k:c:a:d:r:t:h:nfx option
do
  case $option in
    i) DATADIR="$OPTARG"                                                            ;;
    b) BASEFILE="$OPTARG"                                                           ;;
    s) SKETCH=$OPTARG                                                               ;;
    q) Q="$($GAWK -v x=$OPTARG 'BEGIN{printf "%.20f", x+0}' | sed 's/0*$//g')"      ;;
    k) K=$OPTARG                                                                    ;;
    c) CUTOFF="$($GAWK -v x=$OPTARG 'BEGIN{printf "%.20f", x+0}' | sed 's/0*$//g')" ;;
    a) ALPHA="$($GAWK -v x=$OPTARG 'BEGIN{printf "%.20f", x+0}' | sed 's/0*$//g')"  ;;
    f) NFQ=4                                                                        ;;
    n) INFERTREE=false                                                              ;;
    x) BRANCH_SUPPORT=false                                                         ;;
    r) RATCHET=$OPTARG                                                              ;;
    h) CHUNK=$OPTARG                                                                ;;
    t) NPROC=$OPTARG                                                                ;;
    :) echo "missing argument ($OPTARG)" >&2 ; exit 1                               ;;
   \?) echo "invalid option ($OPTARG)"   >&2 ; exit 1                               ;;
  esac
done
if [ "$DATADIR" == "N.O_D.I.R" ];             then echo "genome directory is not specified (option -i)"       >&2 ; exit 1 ; fi
if [ ! -e "$DATADIR" ];                       then echo "genome directory does not exist (option -i)"         >&2 ; exit 1 ; fi
if [ ! -d "$DATADIR" ];                       then echo "$DATADIR is not a directory (option -i)"             >&2 ; exit 1 ; fi
if [ "$BASEFILE" == "N.O_B.A.S.E.F.I.L.E" ];  then echo "basename is not specified (option -b)"               >&2 ; exit 1 ; fi
if [ $(echo "$SKETCH<=0" | bc) -ne 0 ];       then echo "incorrect value: $SKETCH (option -s)"                >&2 ; exit 1 ; fi
if ! [[ $K =~ ^[0-9]+$ ]];                    then echo "incorrect k-mer size: $K (option -k)"                >&2 ; exit 1 ; fi
if ! [[ $Q =~ ^0\.[0-9]+$ ]];                 then echo "incorrect probability q: $Q (option -q)"             >&2 ; exit 1 ; fi
if ! [[ $CUTOFF =~ ^0\.[0-9]+$ ]];            then echo "incorrect cutoff: $CUTOFF (option -c)"               >&2 ; exit 1 ; fi
if ! [[ $ALPHA =~ ^[0-9]+\.[0-9]+$ ]];        then echo "incorrect parameter: $ALPHA (option -a)"             >&2 ; exit 1 ; fi
if ! [[ $RATCHET =~ ^[0-9]+$ ]];              then echo "incorrect ratchet step number: $RATCHET (option -r)" >&2 ; exit 1 ; fi
if ! [[ $NPROC =~ ^[0-9]+$ ]];                then echo "incorrect value: $NPROC (option -t)"                 >&2 ; exit 1 ; fi

### verifying the number of threads #########################################################################
[ $NPROC -le 0 ] && NPROC=2;
echo "$NPROC thread(s)" ;
WAITIME=$($GAWK -v x=$NPROC 'BEGIN{print 1/sqrt(x)}');
						     
### gathering the genome list ###############################################################################
GLIST=$(ls $DATADIR/*.fna $DATADIR/*.fas $DATADIR/*.fa $DATADIR/*.fasta 2>/dev/null);
if   [ $(grep -c -F " " <<<"$GLIST") -ne 0 ]
then
  echo "found blank space(s) within the following file name(s):"                              >&2 ;
  ls $DATADIR/*.fna $DATADIR/*.fas $DATADIR/*.fa $DATADIR/*.fasta 2>/dev/null | grep " "     1>&2 ;
  echo "please avoid any blank space within every file name"                                  >&2 ; exit 1 ;
elif [ $(grep -c -F "'" <<<"$GLIST") -ne 0 ]
then
  echo "found single quote(s) within the following file name(s):"                             >&2 ;
  ls $DATADIR/*.fna $DATADIR/*.fas $DATADIR/*.fa $DATADIR/*.fasta 2>/dev/null | grep -F "'"  1>&2 ;
  echo "please avoid any single quote within every file name"                                 >&2 ; exit 1 ;
elif [ $(grep -c -F "\"" <<<"$GLIST") -ne 0 ]
then
  echo "found double quote(s) within the following file name(s):"                             >&2 ;
  ls $DATADIR/*.fna $DATADIR/*.fas $DATADIR/*.fa $DATADIR/*.fasta 2>/dev/null | grep -F "\"" 1>&2 ;
  echo "please avoid any double quote within every file name"                                 >&2 ; exit 1 ;
elif [ $(grep -c -F ";" <<<"$GLIST") -ne 0 ]
then
  echo "found semicolon(s) within the following file name(s):"                                >&2 ;
  ls $DATADIR/*.fna $DATADIR/*.fas $DATADIR/*.fa $DATADIR/*.fasta 2>/dev/null | grep -F ";"  1>&2 ;
  echo "please avoid any semicolon within every file name"                                    >&2 ; exit 1 ;
elif [ $(grep -c -F "," <<<"$GLIST") -ne 0 ]
then
  echo "found comma(s) within the following file name(s):"                                    >&2 ;
  ls $DATADIR/*.fna $DATADIR/*.fas $DATADIR/*.fa $DATADIR/*.fasta 2>/dev/null | grep -F ","  1>&2 ;
  echo "please avoid any comma within every file name"                                        >&2 ; exit 1 ;
fi
GLIST=$(ls $DATADIR/*.fna $DATADIR/*.fas $DATADIR/*.fa $DATADIR/*.fasta 2>/dev/null | sort);
n=$(echo $GLIST | $GAWK '{print NF}');
if [ $n -lt 4 ]
then
  echo "directory $DATADIR should contain at least 4 files *.fna, *.fas, *.fasta or *.fa"     >&2 ; exit 1 ;
fi
echo "$n taxa" ;

### creating output file names ##############################################################################
ACGT=$BASEFILE.acgt;    # ACGT content of each input genome
OEPL=$BASEFILE.oepl;    # p-distance estimates in OEPL (One Entry Per Line) format 
DMAT=$BASEFILE.d;       # evolutionary distances in PHYLIP square format


#############################################################################################################
#############################################################################################################
#### PREPROCESSING GENOMES                                                                               ####
#############################################################################################################
#############################################################################################################

function ctrl_c() {
  echo -n " process interrupted: deleting files ... " ;
  sleep 5 ;
  for f in $GLIST
  do
    rm -f ${f%.*}.msh ;
  done
  rm -f $ACGT $ACGT.tmp $OEPL $DMAT ;
  echo "[ok]" ;
  exit 1 ;
}
trap  ctrl_c  INT ; 

### estimating ACGT content #################################################################################
rm -f $ACGT ; 
for f in $GLIST
do
  if [ ! -e $f ]
  then
    sleep 1 ; echo "ERROR!" >&2 ;
    sleep 1 ; echo "blank space in  file name starting by:  $f" >&2 ;
    sleep 1 ; echo "please avoid any blank space within every file name" >&2 ;
    sleep 1 ; for f in $GLIST ; do if [ -e $f.nfh ]; then rm -f $f.nfh ; fi ; done
    rm -f $ACGT ;
    exit 1 ;
  fi
  echo "parsing $(basename ${f%.*})" >&2 ;
  $EXEC "grep -v '^>' $f > $f.nfh; a=\$(tr -cd Aa <$f.nfh | wc -c); c=\$(tr -cd Cc <$f.nfh | wc -c); g=\$(tr -cd Gg <$f.nfh | wc -c); t=\$(tr -cd Tt <$f.nfh | wc -c); flock -x $ACGT echo $(basename $f) \$a \$c \$g \$t >> $ACGT; rm -f $f.nfh;" &
  while [ $(jobs -r | wc -l) -ge $NPROC ]; do sleep $WAITIME ; done
done

wait ;

sort $ACGT > $ACGT.tmp ;
mv $ACGT.tmp $ACGT ;

### estimating k-mer and sketch size ########################################################################
[ $K -le 0 ] && K=$($GAWK -v q=$Q '{n=$2+$3+$4+$5;kc=int(log(n*(1-q)/q)/log(4)+0.5);k=(kc>k)?kc:k}END{print k}' $ACGT) && [ $K -le 0 ] && k=19;
echo "k-mer size: $K (q=$Q)" ;
[ $(echo "$SKETCH<=1" | bc) -ne 0 ] && SKETCH=$($GAWK -v s=$SKETCH '{n+=$2+$3+$4+$5}END{printf("%d\n", 1000*int((s*n/NR)/1000))}' $ACGT) && [ $SKETCH -eq 0 ] && SKETCH=10000;
[ $(echo "$SKETCH>1"  | bc) -ne 0 ] && SKETCH=$(echo "($SKETCH+0.5)/1" | bc);
echo "sketch size: $SKETCH" ;

### sketching genomes #######################################################################################
TLIST="" ;
for f in $GLIST
do
  TLIST="$TLIST $(basename ${f%.*})" ;
  echo "sketching $(basename ${f%.*})" >&2 ;
  $EXEC "$MASH sketch -o ${f%.*} -s $SKETCH -k $K $f" &> /dev/null &
  while [ $(jobs -r | wc -l) -ge $NPROC ]; do sleep $WAITIME ; done
done


wait ; 


#############################################################################################################
#############################################################################################################
#### P-DISTANCE ESTIMATES                                                                                ####
#############################################################################################################
#############################################################################################################

### estimating and writing pairwise p-distances #############################################################
echo $TLIST > $OEPL ;
a=($(ls $DATADIR/*.msh | sort));
i=${#a[@]};
while [ $((j=--i)) -ge 0 ]
do
  mi=${a[$i]};
  ti=$(basename ${mi%.*});
  while [ $((--j)) -ge 0 ]
  do
    if [ $CHUNK -eq 1 ] || [ $j -eq 0 ]
    then
      mj=${a[$j]};
      tj=$(basename ${mj%.*});
      echo "estimating p-distance between $ti ($(( $i + 1 ))) and $tj ($(( $j + 1 )))" >&2 ;
      $EXEC "d=\$(timeout 5 $MASH dist -s $SKETCH $mi $mj | $GAWK '{printf(\"%.8f\\n\",1-exp(-\$3))}'); [ -n \"\$d\" ] && flock -x $OEPL echo \"$(( $i + 1 )) $(( $j + 1 )) \$d\" >> $OEPL ;" &
    else
      chk=$CHUNK;
      let j++;
      if [ $j -le $(( $CHUNK - 1 )) ]; then chk=$j; fi
      mlist="";
      c=$chk;
      while [ $((--c)) -ge 0 ]
      do
	let j--;
	mj=${a[$j]};
	mlist="$mlist $mj";
        tj=$(basename ${mj%.*});
        echo "estimating p-distance between $ti ($(( $i + 1 ))) and $tj ($(( $j + 1 )))" >&2 ;
      done
      $EXEC "i=\$(( $i + 1 )); j=\$(( $j + $chk + 1 )); out=; for mj in $mlist; do d=\$(timeout 5 $MASH dist -s $SKETCH $mi \$mj | $GAWK '{printf(\"%.8f\\n\",(\$3==1)?1:1-exp(-\$3))}'); j=\$(( \$j - 1 )); [ -n \"\$d\" ] && out=\"\$out\"\"\$i \$j \$d\\n\"; done ; [ -n \"\$out\" ] && flock -x $OEPL echo -n -e \$out >> $OEPL ;" &
    fi
    while [ $(jobs -r | wc -l) -ge $NPROC ]; do sleep $WAITIME ; done
  done
done

### waiting up to one minute ################################################################################
### if any job is still running (very rare), then killing everything and next re-estimating one by one 
for t in {1..121}
do
  echo -n "." >&2 ;
  [ $(jobs -r | wc -l) -eq 0 ] && break ;
  [ $t -eq 121 ] && jobs -p | xargs kill &>/dev/null;
  sleep 0.5 ;
done

### verifying every p-distance estimate #####################################################################
$GAWK '(NR==1){n=NF;while((j=++i)<=n)while(--j>0)d[i][j]=d[j][i]=-1;next}
              {d[$1][$2]=(d[$2][$1]=$3)}
       END    {i=0;while((j=++i)<=n)while(--j>0)if(d[i][j]<0||d[j][i]<0)print i"\t"j}' $OEPL |
  while read -r i j
  do
    let i--;
    mi=${a[$i]};
    ti=$(basename ${mi%.*});
    let j--;
    mj=${a[$j]};
    tj=$(basename ${mj%.*});
    echo "re-estimating p-distance between $ti ($(( $i + 1 ))) and $tj ($(( $j + 1 )))" >&2 ;
    d=$($MASH dist -s $SKETCH $mi $mj | $GAWK '{printf("%.8f\n",($3==1)?1:1-exp(-$3))}');
    echo "$(( $i + 1 )) $(( $j + 1 )) $d" ;
  done >> $OEPL 

echo " [ok]" >&2 ;


#############################################################################################################
#############################################################################################################
#### EVOLUTIONARY DISTANCE ESTIMATES                                                                     ####
#############################################################################################################
#############################################################################################################

### transforming (if required) p-distances and writing in PHYLIP square format ##############################
### F81/EI transformations can be forced with option -c 0
### to have PC distances (option -a 0) or only p-distances (option -a > 0), F81/EI transformations should be prevented with option -c 1

if [ -n "$($GAWK -v c=$CUTOFF '(NR==1){next}($3>c){print;exit}' $OEPL)" ]
then
  if [ $(echo "$ALPHA<=0" | bc) -ne 0 ]
  then
    ### F81/EI transformation without gamma shape parameter (option -a 0) ###################################
    $GAWK -v p=8 -v f=$NFQ 'function s(x){return x*x}
                            (ARGIND==1)        {++x;sx=$2+$3+$4+$5;if(f==2){a[x]=t[x]=0.5*($2+$5)/sx;c[x]=g[x]=0.5*($3+$4)/sx}else{a[x]=$2/sx;c[x]=$3/sx;g[x]=$4/sx;t[x]=$5/sx}}
                            (ARGIND==2&&FNR==1){while(++n<=NF){m=(m>(l=length(lbl[n]=$n)))?m:l;d[n][n]=0}--n;  print(b=" ")n;x=0.5;while((x*=2)<m)b=b""b;  next}
                            (ARGIND==2)        {d[$1][$2]=(d[$2][$1]=$3)}
                            END                {while(++i<=n){
                                                  printf substr(lbl[i]b,1,m);ai=a[i];ci=c[i];gi=g[i];ti=t[i];j=0;
                                                  while(++j<=n)
                                                    printf(" %."p"f",((dij=d[i][j])==0)?0:((x=1-dij/(1-ai*a[j]-ci*c[j]-gi*g[j]-ti*t[j]))>0)?((s(ai+a[j])+s(ci+c[j])+s(gi+g[j])+s(ti+t[j]))/4-1)*log(x):9.99999999);
                                                  print""}  }' $ACGT $OEPL > $DMAT ;
    echo "F81/EI distances written into $DMAT" ;
  else
    ### F81/EI transformation with gamma shape parameter (option -a > 0) ####################################
    $GAWK -v p=8 -v f=$NFQ -v alp=$ALPHA 'function s(x){return x*x}
                                          (ARGIND==1)        {++x;sx=$2+$3+$4+$5;if(f==2){a[x]=t[x]=0.5*($2+$5)/sx;c[x]=g[x]=0.5*($3+$4)/sx}else{a[x]=$2/sx;c[x]=$3/sx;g[x]=$4/sx;t[x]=$5/sx}}
                                          (ARGIND==2&&FNR==1){while(++n<=NF){m=(m>(l=length(lbl[n]=$n)))?m:l;d[n][n]=0}--n;  print(b=" ")n;x=0.5;while((x*=2)<m)b=b""b;  next}
                                          (ARGIND==2)        {d[$1][$2]=(d[$2][$1]=$3)}
                                          END                {while(++i<=n){
                                                                printf substr(lbl[i]b,1,m);ai=a[i];ci=c[i];gi=g[i];ti=t[i];j=0;
                                                                while(++j<=n)
                                                                  printf(" %."p"f",((dij=d[i][j])==0)?0:((x=1-dij/(1-ai*a[j]-ci*c[j]-gi*g[j]-ti*t[j]))>0)?alp*(1-(s(ai+a[j])+s(ci+c[j])+s(gi+g[j])+s(ti+t[j]))/4)*((x**(-1/alp))-1):9.99999999);
                                                                print""}  }' $ACGT $OEPL > $DMAT ;
    echo "F81/EI distances (a=$ALPHA) written into $DMAT" ;
  fi    
else
  if [ $(echo "$ALPHA<=0" | bc) -ne 0 ]
  then
    ### PC transformation (option -a 0) #####################################################################
    $GAWK -v p=8 '(NR==1){while(++n<=NF){m=(m>(l=length(lbl[n]=$n)))?m:l;d[n][n]=0}--n;  print(b=" ")n;x=0.5;while((x*=2)<m)b=b""b;  next}  
                         {d[$1][$2]=(d[$2][$1]=$3)}
                  END    {while(++i<=n){
                            printf substr(lbl[i]b,1,m);j=0;
                            while(++j<=n)
                              printf(" %."p"f",((dij=d[i][j])==0)?0:(dij<1)?-log(1-dij):9.99999999);
                            print""}  }' $OEPL > $DMAT ;
    echo "PC distances written into $DMAT" ;
  else
    ### p-distance without transformation (option -a > 0) ###################################################
    $GAWK -v p=8 '(NR==1){while(++n<=NF){m=(m>(l=length(lbl[n]=$n)))?m:l;d[n][n]=0}--n;  print(b=" ")n;x=0.5;while((x*=2)<m)b=b""b;  next}  
                         {d[$1][$2]=(d[$2][$1]=$3)}
                  END    {while(++i<=n){
                            printf substr(lbl[i]b,1,m);j=0;
                            while(++j<=n)
                              printf(" %."p"f",((dij=d[i][j])==0)?0:(dij<1)?dij:9.99999999);
                            print""}  }' $OEPL > $DMAT ;
    echo "p-distances written into $DMAT" ;
  fi
fi

### estimating (if required) missing distances, i.e. every missing dij is replaced by min[dix+djx] ##########
miss=$(grep -o -F " 9.99999999" $DMAT | wc -l);
if [ $miss -ne 0 ]
then
  echo "estimating $(( $miss / 2 )) missing distances inside $DMAT" ;
  $GAWK '(NR>1){m=(m>(l=length(lbl[++n]=$(c=j=1))))?m:l;--j;while(++c<=NF){d[n][++j]=$c;(max<$c)&&max=$c}}
         END{while(++i<=n){j=0;
               while(++j<i){x=0;min=2*max;
                 if(d[i][j]!=9.99999999)continue;
                 while(++x<=n)x!=i&&x!=j&&d[i][x]!=9.99999999&&d[j][x]!=9.99999999&&d[i][x]>=0&&d[j][x]>=0&&min>(y=d[i][x]+d[j][x])&&min=y;
                 print "estimated missing distance between "lbl[i]" ("i") and "lbl[j]" ("j"): "min > "/dev/stderr";
                 d[i][j]=d[j][i]=-min}}
             print(b=" ")n;x=0.5;while((x*=2)<m)b=b""b;i=0;
             while(++i<=n){printf substr(lbl[i]b,1,m);j=0;while(++j<=n)printf(" %.8f",((d[i][j]>=0)?d[i][j]:-d[i][j]));print""}  }' $DMAT > $DMAT.tmp ;
  mv $DMAT.tmp $DMAT ;
fi

### deleting all *.msh files ################################################################################
for f in $DATADIR/*.msh ; do rm -f $f ; done


if ! $INFERTREE ; then exit 0 ; fi


#############################################################################################################
#############################################################################################################
#### BME TREE INFERENCE                                                                                  ####
#############################################################################################################
#############################################################################################################

# trap  INT ;
function ctrl_c() {
  echo -n " process interrupted: deleting files ... " ;
  sleep 5 ;
  for x in $(seq 1 $NPROC)
  do
    rm -f $BASEFILE.dd.$x.c $BASEFILE.dd.$x.c_fastme_stat.txt $BASEFILE.dd.$x.n $BASEFILE.dd.$x.n_fastme_stat.txt $BASEFILE.tt.$x.c $BASEFILE.tt.$x.n $BASEFILE.tt.$x.m ;
  done
  rm -f $ACGT $OEPL $BASEFILE.tax $BASEFILE.d $BASEFILE.dd $BASEFILE.dd_fastme_stat.txt $BASEFILE.nwk $BASEFILE.tt ;
  echo "[ok]" ;
  exit 1 ;
}
trap  ctrl_c  INT ; 

echo "searching for the BME phylogenetic tree..." ;

TAXFILE=$BASEFILE.tax;
sed 1d $DMAT | $GAWK '{print "s/@"(++i)"@/"$1"/"}' > $TAXFILE ;

$GAWK '(NR==1){print;next}{printf"@"(++i)"@";j=1;while(++j<=NF)printf" "$j;print""}' $DMAT > $BASEFILE.dd ;
DMAT=$BASEFILE.dd;

BMETREE=$BASEFILE.nwk;  # BME phylogenetic tree in NEWICK format
OUTTREE=$BASEFILE.tt;

### first BME tree inference ################################################################################
$FASTME -i $DMAT -o $OUTTREE -s -f 12 -T 1 &> /dev/null ;
tblo=$(grep -B1 "Performed" $BASEFILE.dd_fastme_stat.txt | sed -n 1p | sed 's/.* //g' | sed 's/\.$//g');
[ -z "$tblo" ] && tblo=$(grep -o ":[0-9\.-]*" $OUTTREE | tr -d :- | paste -sd+ | bc -l | sed 's/^\./0./');
if [ -z "$tblo" ]
then
  tblo=999999;
  echo "  step 0   NaN" >&2 ;
else
  echo "  step 0   $tblo" >&2 ;
  echo "step 0   tbl=$tblo" ;
fi
cp $OUTTREE $BMETREE;
sed -f $TAXFILE $BMETREE > $BMETREE.tmp ; mv $BMETREE.tmp $BMETREE ; # <=> sed -f $TAXFILE -i $BMETREE ;
rm -f $BASEFILE.dd_fastme_stat.txt ;

### ratchet-based search of the BME tree ####################################################################
if [ $RATCHET -gt 0 ]
then
  for x in $(seq 1 $NPROC); do cp $DMAT $DMAT.$x.c ; done

  step=0; eps=0;
  while [ $step -le $RATCHET ]
  do
    step_prev=$step;
    eps_prev=$eps;
    for x in $(seq 1 $NPROC)
    do
      if [ $((++step)) -gt $RATCHET ]; then break; fi

      ### noising evolutionary distances ####################################################################
      v=0.$((++eps))$step; [ $(echo "$v>=0.7" | bc) -eq 1 ] && v=0$(echo "scale=4;$v*$v/1.4" | bc -l);
      $GAWK -v v=$v -v s=$step 'BEGIN  {srand(s)}
                                (NR==1){n=$0;next}  {lbl[++i]=$1;d[i][i]=0;j=0;f=1;while(++f<=i){++j;d[i][j]=(d[j][i]=($f*(1-v)+2*v*$f*rand()))}}
                                END    {print" "n;i=0;while(++i<=n){printf lbl[i];j=0;while(++j<=n){printf(" %.8f",d[i][j])}print""}}' $DMAT.$x.c > $DMAT.$x.n ;

      ### ratchet-search tree search ########################################################################
      $EXEC "$FASTME -i $DMAT.$x.n -u $OUTTREE -o $OUTTREE.$x.n -nB -s -T 1 ; sed 's/:-/:/g' $OUTTREE.$x.n > $OUTTREE.$x.m ; $FASTME -i $DMAT.$x.c -u $OUTTREE.$x.m -o $OUTTREE.$x.c -s -T 1 -f 12 ; rm -f $DMAT.$x.n $DMAT.$x.m $DMAT.$x.n_fastme_stat.txt $OUTTREE.$x.n $OUTTREE.$x.m ;" &> /dev/null &
    done

    while [ $(jobs -r | wc -l) -gt 0 ]; do sleep $WAITIME ; done

    for x in $(seq 1 $NPROC)
    do
      if [ $((++step_prev)) -gt $RATCHET ]; then break; fi
      v=0.$((++eps_prev))$step_prev; [ $(echo "$v>=0.7" | bc) -eq 1 ] && v=0$(echo "scale=4;$v*$v/1.44" | bc -l);
    
      tbl=$(grep -B1 "Performed" $DMAT.$x.c_fastme_stat.txt | sed -n 1p | sed 's/.* //g' | sed 's/\.$//g');
      rm -f $DMAT.$x.c_fastme_stat.txt ;
      out=" ";
      [ -z "$tbl" ] && tbl=$(grep -o ":[0-9\.-]*" $OUTTREE.$x.c | tr -d :- | paste -sd+ | bc | sed 's/^\./0./') && out="+";
      [ -z "$tbl" ] && tbl="NaN";
      echo -n "$out step $step_prev   $tbl" >&2 ; 
      if [ "$tbl" == "NaN" ] || [ $(echo "$tbl<$tblo" | bc) -eq 0 ]
      then
        rm -f $OUTTREE.$x.c ;
        echo "  (epsilon=$v)" >&2 ;
      else
        tblo=$tbl;
        mv $OUTTREE.$x.c $OUTTREE ;
        cp $OUTTREE $BMETREE;
        sed -f $TAXFILE $BMETREE > $BMETREE.tmp ; mv $BMETREE.tmp $BMETREE ; # <=> sed -f $TAXFILE -i $BMETREE ;
        echo " *  (epsilon=$v)" >&2 ;
        echo "step $step_prev   tbl=$tbl";
        eps=0;
      fi
    done
  done

  for x in $(seq 1 $NPROC); do rm -f $DMAT.$x.c ; done
fi

#############################################################################################################
#############################################################################################################
#### REQ CONFIDENCE VALUES AT BRANCHES                                                                   ####
#############################################################################################################
#############################################################################################################

if $BRANCH_SUPPORT
then
  echo -n "estimating branch supports ... " >&2 ;
  $REQ $BASEFILE.d $BMETREE $OUTTREE ;
  echo "[ok]" >&2 ;
  mv $OUTTREE $BMETREE ;
  echo "BME tree (tbl=$tblo) with branch supports written into $BMETREE" ;
else
  echo "BME tree (tbl=$tblo) written into $BMETREE" ;
fi

rm -f $DMAT $TAXFILE $OUTTREE ;

exit ;

