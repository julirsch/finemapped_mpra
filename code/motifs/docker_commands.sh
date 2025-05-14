sudo docker build -t jcu_motifliq1 .
sudo docker run -i -t --cpus=220 -v `pwd`:/data a9bd6b4cf894 /bin/bash 

for i in {{1..22},X}; do
echo ${i}
/usr/local/bin/motif_liquidator /data/hocomoco.meme /data/mpra_chr${i}.fasta -o mpra_chr${i}_hocomoco.out
gzip mpra_chr${i}_hocomoco.out &
done

for i in {{1..22},X}; do
echo ${i}
/usr/local/bin/motif_liquidator /data/jaspar.meme /data/mpra_chr${i}.fasta -o mpra_chr${i}_jaspar.out
gzip mpra_chr${i}_jaspar.out &
done

for i in {{1..22},X}; do
echo ${i}
/usr/local/bin/motif_liquidator /data/vierstra.meme /data/mpra_chr${i}.fasta -o mpra_chr${i}_vierstra.out
gzip mpra_chr${i}_vierstra.out &
done

for i in {{1..22},X}; do
echo ${i}
/usr/local/bin/motif_liquidator /data/cisbp.meme /data/mpra_chr${i}.fasta -o mpra_chr${i}_cisbp.out
gzip mpra_chr${i}_cisbp.out &
done

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

sudo docker system prune -a -f

/usr/local/bin/motif_liquidator /data/vierstra.meme /data/mpra_chr${i}.fasta -o mpra_chr${i}_vierstra.out
gzip mpra_chr${i}_vierstra.out &
/usr/local/bin/motif_liquidator /data/cisbp.meme /data/mpra_chr${i}.fasta -o mpra_chr${i}_cisbp.out
gzip mpra_chr${i}_cisbp.out &
/usr/local/bin/motif_liquidator /data/jaspar.meme /data/mpra_chr${i}.fasta -o mpra_chr${i}_jaspar.out
gzip mpra_chr${i}_jaspar.out &
