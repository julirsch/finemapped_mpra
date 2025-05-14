#!/bin/bash

for i in {{1..22},X}; do
sudo docker cp 30f14884f0c7:/mpra_chr${i}_hocomoco.out.gz /mnt/sdb/gtex_mpra/data/fimo/.
done

for i in {{1..22},X}; do
sudo docker cp 30f14884f0c7:/mpra_chr${i}_jaspar.out.gz /mnt/sdb/gtex_mpra/data/fimo/.
done

for i in {{1..22},X}; do
sudo docker cp 30f14884f0c7:/mpra_chr${i}_vierstra.out.gz /mnt/sdb/gtex_mpra/data/fimo/.
done

for i in {{1..22},X}; do
sudo docker cp 30f14884f0c7:/mpra_chr${i}_cisbp.out.gz /mnt/sdb/gtex_mpra/data/fimo/.
done

