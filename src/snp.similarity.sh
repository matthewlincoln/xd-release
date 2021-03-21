#!/bin/bash

# Get input and output directories from the command line:
input_direc=$1
ped_direc=$2
output_direc=$3

mkdir -p ${output_direc}/ced \
         ${output_direc}/ibd \
         ${output_direc}/ms \
         ${output_direc}/sle_g \
         ${output_direc}/sle_o \
         ${output_direc}/t1d \
         ${output_direc}/t1d_asp

# Calculate similarity statistics for SNPs mapping to the same location in the CeD cohort:

# Divide duplicate file into smaller blocks of 100 SNPs each:
sed -n '1,100 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.1
sed -n '101,200 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.2
sed -n '201,300 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.3
sed -n '301,400 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.4
sed -n '401,500 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.5
sed -n '501,600 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.6
sed -n '601,700 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.7
sed -n '701,800 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.8
sed -n '801,900 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.9
sed -n '901,1000 p' ${input_direc}/ced.dup.snps.txt > ${output_direc}/ced/ced.dup.snps.txt.10

for i in `seq 1 10`;
do
  sbatch -J ced.dups.${i} \
         -o ${output_direc}/ced/ced.dup.similarity.${i}.out \
         -e ${output_direc}/ced/ced.dup.similarity.${i}.err \
         src/dup.similarity.sh ${output_direc}/ced/ced.dup.snps.txt.${i} \
         ${ped_direc}/ced.original.ped \
         ${ped_direc}/ced.original.map \
         ${output_direc}/ced
done


# Calculate similarity statistics for SNPs mapping to the same location in the IBD cohort:

# Divide duplicate file into smaller blocks of 100 SNPs each:
sed -n '1,100 p' ${input_direc}/ibd.dup.snps.txt > ${output_direc}/ibd/ibd.dup.snps.txt.1
sed -n '101,200 p' ${input_direc}/ibd.dup.snps.txt > ${output_direc}/ibd/ibd.dup.snps.txt.2
sed -n '201,300 p' ${input_direc}/ibd.dup.snps.txt > ${output_direc}/ibd/ibd.dup.snps.txt.3
sed -n '301,400 p' ${input_direc}/ibd.dup.snps.txt > ${output_direc}/ibd/ibd.dup.snps.txt.4
sed -n '401,500 p' ${input_direc}/ibd.dup.snps.txt > ${output_direc}/ibd/ibd.dup.snps.txt.5
sed -n '501,600 p' ${input_direc}/ibd.dup.snps.txt > ${output_direc}/ibd/ibd.dup.snps.txt.6
sed -n '601,700 p' ${input_direc}/ibd.dup.snps.txt > ${output_direc}/ibd/ibd.dup.snps.txt.7
sed -n '701,800 p' ${input_direc}/ibd.dup.snps.txt > ${output_direc}/ibd/ibd.dup.snps.txt.8
sed -n '801,900 p' ${input_direc}/ibd.dup.snps.txt > ${output_direc}/ibd/ibd.dup.snps.txt.9

for i in `seq 1 9`;
do
  sbatch -J ibd.dups.${i} \
         -o ${output_direc}/ibd/ibd.dup.similarity.${i}.out \
         -e ${output_direc}/ibd/ibd.dup.similarity.${i}.err \
         src/dup.similarity.sh ${output_direc}/ibd/ibd.dup.snps.txt.${i} \
         ${ped_direc}/ibd.original.ped \
         ${ped_direc}/ibd.original.map \
         ${output_direc}/ibd
done


# Calculate similarity statistics for SNPs mapping to the same location in the MS cohort:

# Since there are so few duplicates in the MS cohort, there is no need to divide them up
sbatch -J ms.dups \
       -o ${output_direc}/ms/ms.dup.similarity.out \
       -e ${output_direc}/ms/ms.dup.similarity.err \
       src/dup.similarity.sh ${input_direc}/ms.dup.snps.txt \
       ${ped_direc}/ms.original.ped \
       ${ped_direc}/ms.original.map \
       ${output_direc}/ms


# Calculate similarity statistics for SNPs mapping to the same location in the SLE_G cohort:

# Divide duplicate file into smaller blocks of 100 SNPs each:
sed -n '1,100 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.1
sed -n '101,200 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.2
sed -n '201,300 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.3
sed -n '301,400 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.4
sed -n '401,500 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.5
sed -n '501,600 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.6
sed -n '601,700 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.7
sed -n '701,800 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.8
sed -n '801,900 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.9
sed -n '901,1000 p' ${input_direc}/sle_g.dup.snps.txt > ${output_direc}/sle_g/sle_g.dup.snps.txt.10

for i in `seq 1 10`;
do
  sbatch -J sleg.dups.${i} \
         -o ${output_direc}/sle_g/sle_g.dup.similarity.${i}.out \
         -e ${output_direc}/sle_g/sle_g.dup.similarity.${i}.err \
         src/dup.similarity.sh ${output_direc}/sle_g/sle_g.dup.snps.txt.${i} \
         ${ped_direc}/sle_g.original.ped \
         ${ped_direc}/sle_g.original.map \
         ${output_direc}/sle_g
done


# Calculate similarity statistics for SNPs mapping to the same location in the SLE_O cohort:

# Divide duplicate file into smaller blocks of 100 SNPs each:
sed -n '1,100 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.1
sed -n '101,200 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.2
sed -n '201,300 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.3
sed -n '301,400 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.4
sed -n '401,500 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.5
sed -n '501,600 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.6
sed -n '601,700 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.7
sed -n '701,800 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.8
sed -n '801,900 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.9
sed -n '901,1000 p' ${input_direc}/sle_o.dup.snps.txt > ${output_direc}/sle_o/sle_o.dup.snps.txt.10

for i in `seq 1 10`;
do
  sbatch -J sleo.dups.${i} \
         -o ${output_direc}/sle_o/sle_o.dup.similarity.${i}.out \
         -e ${output_direc}/sle_o/sle_o.dup.similarity.${i}.err \
         src/dup.similarity.sh ${output_direc}/sle_o/sle_o.dup.snps.txt.${i} \
         ${ped_direc}/sle_o.original.ped \
         ${ped_direc}/sle_o.original.map \
         ${output_direc}/sle_o
done


# Calculate similarity statistics for SNPs mapping to the same location in the T1D cohort:

# Divide duplicate file into smaller blocks of 100 SNPs each:
sed -n '1,100 p' ${input_direc}/t1d.dup.snps.txt > ${output_direc}/t1d/t1d.dup.snps.txt.1
sed -n '101,200 p' ${input_direc}/t1d.dup.snps.txt > ${output_direc}/t1d/t1d.dup.snps.txt.2
sed -n '201,300 p' ${input_direc}/t1d.dup.snps.txt > ${output_direc}/t1d/t1d.dup.snps.txt.3
sed -n '301,400 p' ${input_direc}/t1d.dup.snps.txt > ${output_direc}/t1d/t1d.dup.snps.txt.4
sed -n '401,500 p' ${input_direc}/t1d.dup.snps.txt > ${output_direc}/t1d/t1d.dup.snps.txt.5
sed -n '501,600 p' ${input_direc}/t1d.dup.snps.txt > ${output_direc}/t1d/t1d.dup.snps.txt.6
sed -n '601,700 p' ${input_direc}/t1d.dup.snps.txt > ${output_direc}/t1d/t1d.dup.snps.txt.7
sed -n '701,800 p' ${input_direc}/t1d.dup.snps.txt > ${output_direc}/t1d/t1d.dup.snps.txt.8
sed -n '801,900 p' ${input_direc}/t1d.dup.snps.txt > ${output_direc}/t1d/t1d.dup.snps.txt.9

for i in `seq 1 9`;
do
  sbatch -J t1d.dups.${i} \
         -o ${output_direc}/t1d/t1d.dup.similarity.dups.${i}.out \
         -e ${output_direc}/t1d/t1d.dup.similarity.dups.${i}.err \
         src/dup.similarity.sh ${output_direc}/t1d/t1d.dup.snps.txt.${i} \
         ${ped_direc}/t1d.original.ped \
         ${ped_direc}/t1d.original.map \
         ${output_direc}/t1d
done

# Calculate similarity statistics for SNPs mapping to the same location in the T1D_ASP cohort:

# Divide duplicate file into smaller blocks of 100 SNPs each:
sed -n '1,100 p' ${input_direc}/t1d_asp.dup.snps.txt > ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.1
sed -n '101,200 p' ${input_direc}/t1d_asp.dup.snps.txt > ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.2
sed -n '201,300 p' ${input_direc}/t1d_asp.dup.snps.txt > ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.3
sed -n '301,400 p' ${input_direc}/t1d_asp.dup.snps.txt > ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.4
sed -n '401,500 p' ${input_direc}/t1d_asp.dup.snps.txt > ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.5
sed -n '501,600 p' ${input_direc}/t1d_asp.dup.snps.txt > ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.6
sed -n '601,700 p' ${input_direc}/t1d_asp.dup.snps.txt > ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.7
sed -n '701,800 p' ${input_direc}/t1d_asp.dup.snps.txt > ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.8
sed -n '801,900 p' ${input_direc}/t1d_asp.dup.snps.txt > ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.9

for i in `seq 1 9`;
do
  sbatch -J t1d_asp.dups.${i} \
         -o ${output_direc}/t1d_asp/t1d_asp.dup.similarity.dups.${i}.out \
         -e ${output_direc}/t1d_asp/t1d_asp.dup.similarity.dups.${i}.err \
         src/dup.similarity.sh ${output_direc}/t1d_asp/t1d_asp.dup.snps.txt.${i} \
         ${ped_direc}/t1d_asp.original.ped \
         ${ped_direc}/t1d_asp.original.map \
         ${output_direc}/t1d_asp
done
