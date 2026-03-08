#!/bin/bash
#SBATCH --mem 1G
#SBATCH -t 3-00:00:00

#source /home/fs70587/lijuan2/.bashrc

echo "start: `date`"

#wd=/home/fs70587/lijuan2/polygenicAdaptation/
#code=$wd/code 
##Note: sim*.py are just files to import the main simulation functions for easy to read the parameters, not necessary. The main functions could be called directly. 

#:<<A
###### Figure 1; single-locus underdominant model.
#initFrq=0.5
#data=$wd/data/underDomi_replicates
#Ns=$1
#Nm=$2

#for Ns in 0.1 1 5
#do
#	for Nm in 0 0.03 0.05 0.1 0.2 0.3 0.5 0.7 1 2 5 10
#	do
#		for r in $(seq 10)
#		do
#		python3 $code/simUnderDomi.py  --max_generations 500001 --Ns $Ns --Nm $Nm --init_allele_freq $initFrq --outfile $data/frq_$initFrq\_Nm$Nm\_Ns$Ns\_r$r
#		done
#	done
#done
#A

#:<<B
##### Figure 2; stabilizing selection; equal effect sizes;
#data=$wd/data/stabSel_replicates/
#Ns=$1
#Nm=$2
#r=$3

#for Ns in 0.1 1 5
#do
#	for Nm in 0 0.03 0.05 0.1 0.2 0.3 0.5 0.7 1 2 5 10
#	do
#		for r in 1 2 3 4
#		do
#			echo $Nm
#			if [ $(echo $Nm '<0.5' | bc -l) -eq 0 ]
#			then
#				echo "ge"
#				python3 $code/simStabSelV1B1.py --Nu 0.01 --Ns $Ns --Nm $Nm --num_demes 100 --num_loci 200 --deme_size 200   --num_generations 100002  --output_file  $data/B1/het_trait_frqB1_Nu0.01_Nm$Nm\_Ns$Ns\_N200_D100_L200_r$r
#			else
#				echo "le"
#				python3 $code/simStabSelV1B1.py --Nu 0.01 --Ns $Ns --Nm $Nm --num_demes 100 --num_loci 200 --deme_size 200   --num_generations 500002  --output_file  $data/B1/het_trait_frqB1_Nu0.01_Nm$Nm\_Ns$Ns\_N200_D100_L200_r$r
#			fi
#		done
#	done
#done
#B

#:<<C
##### Figure 3; stabilizing selection; different effect sizes. 
###three different effect sizes. (figure 3a,3b,S6)
#data=$wd/data/stabSel_diffEff/
#Nm=$1
#r=$1

#for Nm in 0 0.03 0.05 0.1 0.2 0.3 0.5 0.7 1 2 5 10
#do
#	for r in 1 2 3 4
#	do
#		echo $Nm
#		echo $r
#		numLoc=200
#		Vs=1
#		N=200
#		if [ $(echo $Nm '<0.5' | bc -l) -eq 0 ]
#		then
#			echo "ge"
#			python3 $code/simStabSel_diffEff.py  --subpopulation_size $N --Nm $Nm --Vs $Vs --num_loci $numLoc --num_generations 100002  --output_file $data/3eff/het_trait_frqB1_Nm$Nm\_Vs$Vs\_N$N\_D100_L$numLoc\_r$r
#		else
#			echo "le"
#			python3 $code/simStabSel_diffEff.py  --subpopulation_size $N --Nm $Nm --Vs $Vs --num_loci $numLoc --num_generations 500002  --output_file $data/3eff/het_trait_frqB1_Nm$Nm\_Vs$Vs\_N$N\_D100_L$numLoc\_r$r
#		fi
#	done
#done
#C

#:<<D
###### mixed effect sizes (figure 3c,3d,4e,4f)
#Nm=10
#D=100
#r=$1
#echo "Figure 3c,d $Nm"
#for Nm in 0 0.1 0.5 2
#do
#	for r in 1
#	do
#		if [ $r -eq 1 ]
#		then
#			python3 $code/simStabSelMixEffSize.py --num_generations 500000 --Nm $Nm --Vs 1  --num_loci 200 --deme_size 200 --num_demes $D --Nu 0.01  --output_file $data/mixEff/Nu0.01_Vs1_aExpNuL_Nm$Nm\_L200_N200_D$D\_r$r --prop 1
#		else
#			python3 $code/simStabSelMixEffSize.py --num_generations 500000 --Nm $Nm --Vs 1  --num_loci 200 --deme_size 200 --num_demes $D --Nu 0.01  --output_file $data/mixEff/Nu0.01_Vs1_aExpNuL_Nm$Nm\_L200_N200_D$D\_r$r --prop 1 --effSize_file $data/mixEff/Nu0.01_Vs1_aExpNuL_Nm$Nm\_L200_N200_D$D\_r1.freq
#		fi
#	done
#done
#D


#:<<E
##### Figure 4; equal effect sizes. (figure 4a-d,S3)
#data=$wd/data/stabSel_replicates/Fv
#Ns=5
#Nm=$1
#r=$2

#echo "Ns: $Ns; Nm: $Nm; r: $r"

#for Ns in 0.2 1 5
#do
#	for Nm in 0 0.25 10 0.089 2
#	do
#		for r in 1 2 3 4
#		do
#			if [ $(echo $Nm '<0.5' | bc -l) -eq 0 ]
#			then
#				python3 $code/simStabSelV1Frq.py --Nu 0.01 --Ns $Ns --Nm $Nm --num_demes 100 --num_loci 200 --deme_size 200 --num_generations 100001 --output_file $data/frqB1_Nu0.01_Nm$Nm\_Ns$Ns\_N200_D100_L200_r$r
#			else
#				python3 $code/simStabSelV1Frq.py --Nu 0.01 --Ns $Ns --Nm $Nm --num_demes 100 --num_loci 200 --deme_size 200 --num_generations 500001 --output_file $data/frqB1_Nu0.01_Nm$Nm\_Ns$Ns\_N200_D100_L200_r$r
#			fi

#		done
#	done
#done
#E

#:<<F
#### Figure S2; the effect of number of demes
#data=$wd/data/underDomi_numDeme/
#Nm=0.5
#r=8
#D=25

##Fixed Nu
#for D in 25 50 250 500
#do
#	for r in  $(seq 10)
#	do
#		python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes $D  --pop_size 100 --Ns 0.1 --Nu 0.01 --outfile $data/frq0.5_Nm$Nm\_D$D\_N100_Ns0.1_Nu0.01_r$r
#	done
#done

## Fixed NuD
#for Nm in 0.5 2
#for r in $(seq 10)
#do

#echo "NuD: 5; r: $r"
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 25  --pop_size 2000 --Ns 2 --Nu 0.2 --outfile $data/frq0.5_Nm$Nm\_D25_N2000_Ns2_Nu0.2_r$r
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 50  --pop_size 1000 --Ns 1   --Nu 0.1 --outfile $data/frq0.5_Nm$Nm\_D50_N1000_Ns1_Nu0.1_r$r
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 100 --pop_size 500  --Ns 0.5 --Nu 0.05 --outfile $data/frq0.5_Nm$Nm\_D100_N500_Ns0.5_Nu0.05_r$r
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 250 --pop_size 200  --Ns 0.2 --Nu 0.02 --outfile $data/frq0.5_Nm$Nm\_D250_N200_Ns0.2_Nu0.02_r$r
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 500 --pop_size 100  --Ns 0.1 --Nu 0.01 --outfile $data/frq0.5_Nm$Nm\_D500_N100_Ns0.1_Nu0.01_r$r

#echo "NuD: 1; r: $r"
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 25  --pop_size 400 --Ns 0.4 --Nu 0.04 --outfile $data/frq0.5_Nm$Nm\_D25_N400_Ns0.4_Nu0.04_r$r
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 50  --pop_size 200 --Ns 0.2   --Nu 0.02 --outfile $data/frq0.5_Nm$Nm\_D50_N200_Ns0.2_Nu0.02_r$r
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 100 --pop_size 100  --Ns 0.1 --Nu 0.01 --outfile $data/frq0.5_Nm$Nm\_D100_N100_Ns0.1_Nu0.01_r$r
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 250 --pop_size 40  --Ns 0.04 --Nu 0.004 --outfile $data/frq0.5_Nm$Nm\_D250_N40_Ns0.04_Nu0.004_r$r
#python3 $code/simUnderDomi.py --max_generations 500002 --init_allele_freq 0.5 --Nm $Nm --num_demes 500 --pop_size 20  --Ns 0.02 --Nu 0.002 --outfile $data/frq0.5_Nm$Nm\_D500_N20_Ns0.02_Nu0.002_r$r
#done
#done

#F

#:<<G
#### Figure S4; isolated population; to find out the effect of NuL to genic variance
#data=$wd/data/stabSel_numLoci/
#Nu=0.01
#Ns=5
#numLoc=$1
#r=$2

# echo "Figure S4"
#for numLoc in $(seq 100 100 1000)
#do
#	demeSize=$numLoc
#	for r in $(seq 11 20)
#	do
#		 python3 $code/simStabSelV1.py --num_demes 1  --deme_size $demeSize  --Nu $Nu --num_loci $numLoc --Ns $Ns --Nm 0 --num_generations 50001  --output_file $wd/data/stabSel_numLoci/initBeta/het_trait_Nu$Nu\_Ns$Ns\_N$demeSize\_D1_L$numLoc\_r$r 
#	done
#done

#G

#:<<H
#### Figure S5; the effect of intial condition on equilibration
#data=$wd/data/diffInit/

#for Ns in 1 5
#do
#	for Nm in 0.1 0.2 0.5
#	do
#		for init in 0.5 0
#		do
#			echo $init
#			python $code/simUnderDomi.py --max_generations 1000001 --Ns 5 --Nm $Nm --num_demes 100 --init_allele_freq $init --pop_size 100 --outfile $data/init$init\_frq_Nm$Nm\_Ns$Ns\_d100_N100 &
#	       	done
		
		### same intial condiation for all demes
#			python3 $code/simStabSelV1.py --Nu 0.01 --Ns $Ns --Nm $Nm --num_demes 100 --num_loci 200 --deme_size 200   --num_generations 1000001  --output_file  $data/het_trait_frqB0_Nu0.01_Nm$Nm\_Ns$Ns\_N200_D100_L200 &
		### different intial condition for each deme
#			python3 $code/simStabSelV1B1.py --Nu 0.01 --Ns $Ns --Nm $Nm --num_demes 100 --num_loci 200 --deme_size 200   --num_generations 1000001  --output_file  $data/het_trait_frqB1_Nu0.01_Nm$Nm\_Ns$Ns\_N200_D100_L200 &
#        done
#done
#H



echo "end: `date`"
