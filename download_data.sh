#!/bin/bash
# Downloads necessary data files for the Minimap2 tutorial from OSG Public Storage to the local environment.
mkdir -p /ospool/guest-ap/data/jovyan/tutorial-minimap2/{software,inputs}
mkdir -p ~/tutorial-minimap2/{inputs, outputs}

#osdf object get /ospool/uc-shared/public/osg-training/tutorial-minimap2/humpback_whale_genome_09OCT2025_v1.fasta ./inputs/
#osdf object get /ospool/uc-shared/public/osg-training/tutorial-minimap2/h_whale_sd_0001_sub.fastq.gz ./inputs/
#osdf object get /ospool/uc-shared/public/osg-training/tutorial-minimap2/minimap2_08OCT2025_v1.sif /ospool/guest-ap/data/jovyan/tutorial-minimap2/software/
#osdf object get /ospool/uc-shared/public/osg-training/tutorial-minimap2/humpback_whale_genome_09OCT2025_v1.fasta.mmi /ospool/guest-ap/data/jovyan/tutorial-minimap2/inputs/

osdf object get /ospool/uc-shared/public/osg-training/tutorial-minimap2/minimap2_08OCT2025_v1.sif /ospool/guest-ap/data/jovyan/tutorial-minimap2/software/ > /dev/null
echo "A pre-built container of Minimap2 is at /ospool/guest-ap/data/jovyan/tutorial-minimap2/software/minimap2_08OCT2025_v1.sif ✅"
echo ""

osdf object get /ospool/uc-shared/public/osg-training/tutorial-minimap2/humpback_whale_ref_genome.fasta.mmi /ospool/guest-ap/data/jovyan/tutorial-minimap2/inputs/ > /dev/null
echo "A pre-indexed reference genome is at /ospool/guest-ap/data/jovyan/tutorial-minimap2/inputs/humpback_whale_ref_genome.fasta.mmi ✅"
echo ""

osdf object get /ospool/uc-shared/public/osg-training/tutorial-minimap2/h_whale_sd_0001_sub.fastq.tar.gz?pack=tar.gz ./inputs/ > /dev/null
mv inputs/h_whale_sd_0001_sub.fastq inputs/humpback_whale_reads.fastq
echo "A copy of the WGS reads is at ./inputs/h_whale_sd_0001_sub.fastq ✅"
echo "Data download complete. ✅"