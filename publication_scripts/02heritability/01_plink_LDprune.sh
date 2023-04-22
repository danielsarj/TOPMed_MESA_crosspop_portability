for t in PBMC 
do
	for p in AFA EUR HIS CHN
	do
		for r in 0.5 0.2
		do

        	plink --bfile /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.hg38 --indep-pairwise 500 50 ${r} --ld-xchr 1 --out /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.${r}.hg38 && plink --bfile /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.hg38 --extract /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.${r}.hg38.prune.in --make-bed --out /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.${r}.hg38.pruned && gzip /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.${r}*

		done
	done
done

for t in Mono Tcell
do
	for p in AFA EUR HIS
	do
		for r in 0.5 0.2
		do

        	plink --bfile /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.hg38 --indep-pairwise 500 50 ${r} --ld-xchr 1 --out /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.${r}.hg38 && plink --bfile /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.hg38 --extract /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.${r}.hg38.prune.in --make-bed --out /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.${r}.hg38.pruned && gzip /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.${t}.${p}.${r}*

		done
	done
done

gzip /home/daniel/MESA_heritability/plink_files/* && gunzip /home/daniel/MESA_heritability/plink_files/*.txt.gz

