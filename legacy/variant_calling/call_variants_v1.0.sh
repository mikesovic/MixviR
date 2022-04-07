if [[ ! -d vcfs ]]; then
    mkdir vcfs
fi

if [[ ! -d AD ]]; then
    mkdir AD
fi

if [[ ! -d depths ]]; then
    mkdir depths
fi

if [[ ! -d consensus ]]; then
    mkdir consensus
fi

if [[ ! -d Rout ]]; then
    mkdir Rout
fi

for i in `ls *.bam`
do
qsub call_variants.job ${i}
done
