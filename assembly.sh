#! /bin/bash

# Default values
THREADS=1
RACON_ITER=1
PILON_ITER=1
WTDGB2_PRESET=""
ASSEMBLER=""
TECH=""
GENOME_SIZE=""
OUTDIR=""
LONG_READS=""
READ1=""
READ2=""

if [[ $TECH = "nanopore" ]]; then
MINIMAP2_PRESET="ava-ont"
else
MINIMAP2_PRESET="ava-pb"
fi

usage() {
    echo "Usage: $0 [options]"
    echo "  -l <long reads file>                Long reads FASTQ (gzipped or not)"
    echo "  -1 <short read R1>                  Paired-end short reads R1"
    echo "  -2 <short read R2>                  Paired-end short reads R2"
    echo "  -asm <assembler>                    Assembler: canu or wtdbg2"
    echo "  -tech <technology>                  Technology: pacbio or nanopore"
    echo "  -g|--genome-size <genome size>      Genome size (e.g., 5m, 2.6g)"
    echo "  -o|--outdir <output directory>      Output directory"
    echo "  -t|--threads <threads>              Number of threads (default: 1)"
    echo "  --racon-iter <N>                    Racon polishing iterations (default: 0)"
    echo "  --pilon-iter <N>                    Pilon polishing iterations (default: 0)"
    echo "  --wtdbg2-preset <preset(s)>         Presets for wtdbg2 (comma-separated)"
    # --wtdbg2-preset: Preset for wtdbg2, comma seperable i.e. `ont,preset2` for ont reads and expected genome size <1G`
    #       preset1/rsII/rs: -p 21 -S 4 -s 0.05 -L 5000
    #       preset2: -p 0 -k 15 -AS 2 -s 0.05 -L 5000
    #       preset3: -p 19 -AS 2 -s 0.05 -L 5000
    #       sequel/sq
    #       nanopore/ont:
    #           (genome size < 1G: preset2) -p 0 -k 15 -AS 2 -s 0.05 -L 5000
    #           (genome size >= 1G: preset3) -p 19 -AS 2 -s 0.05 -L 5000
    #       preset4/corrected/ccs: -p 21 -k 0 -AS 4 -K 0.05 -s 0.5
    exit 1
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -l) LONG_READS="$2"; shift ;;
        -1) READ1="$2"; shift ;;
        -2) READ2="$2"; shift ;;
        -asm) ASSEMBLER="$2"; shift ;;
        -tech) TECH="$2"; shift ;;
        -g|--genome-size) GENOME_SIZE="$2"; shift ;;
        -o|--outdir) OUTDIR="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        --racon-iter) RACON_ITER="$2"; shift ;;
        --pilon-iter) PILON_ITER="$2"; shift ;;
        --wtdbg2-preset) WTDGB2_PRESET="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

HALF_THREADS=$(printf "%.0f" "$(echo "$THREADS * 0.5" | bc -l)")

# Validate required arguments
if [[ -z "$LONG_READS" || -z "$ASSEMBLER" || -z "$TECH" || -z "$OUTDIR" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

mkdir -p ${OUTDIR}/01_assembly
contigs_fasta="${OUTDIR}/01_assembly/assembly.contigs.fasta"
case $ASSEMBLER in
    "canu")
        if [[ -z "$GENOME_SIZE" ]]; then
            echo "Error: Canu requires an estimated genome size."
            echo "       -g|--genome-size <genome size>      Genome size (e.g., 5m, 2.6g)"
            exit 1
        fi
        canu \
            -p assembly \
            -d ${OUTDIR}/01_assembly \
            genomeSize=${GENOME_SIZE} \
            -${TECH}-raw ${LONG_READS}
        ;;
    "wtdbg2")
        if [[ -z "$WTDGB2_PRESET" ]]; then
            echo "Error: Please specify --wtdbg2-preset"
            echo "For example: ont,preset2 for ont reads and expected genome size <1G"
            echo "- preset1/rsII/rs: -p 21 -S 4 -s 0.05 -L 5000"
            echo "- preset2: -p 0 -k 15 -AS 2 -s 0.05 -L 5000"
            echo "- preset3: -p 19 -AS 2 -s 0.05 -L 5000"
            echo "- sequel/sq"
            echo "- nanopore/ont:"
            echo -e "\t(genome size < 1G: preset2) -p 0 -k 15 -AS 2 -s 0.05 -L 5000"
            echo -e "\t(genome size >= 1G: preset3) -p 19 -AS 2 -s 0.05 -L 5000"
            echo "- preset4/corrected/ccs: -p 21 -k 0 -AS 4 -K 0.05 -s 0.5"
            exit 1
        fi
        if [[ -n "$GENOME_SIZE" ]]; then
            ADDITIONAL_wtdbg2_PARAMS="-g ${GENOME_SIZE} ${ADDITIONAL_wtdbg2_PARAMS}"
        fi
        wtdbg2 \
            -x ${WTDGB2_PRESET} \
            -i ${LONG_READS} \
            -o ${OUTDIR}/01_assembly/assembly \
            -t ${THREADS} -f ${ADDITIONAL_wtdbg2_PARAMS} 2>&1 | tee ${OUTDIR}/01_assembly/wtdbg2.log
        exit_code=$?
        if [[ "$exit_code" -ne 0 ]]; then
            echo "ERROR assembly with wtdbg2."
            exit $exit_code
        fi 
        wtpoa-cns -i ${OUTDIR}/01_assembly/assembly.ctg.lay.gz -o ${contigs_fasta} -t ${THREADS} -f
        ;;
    "spades")
        if [[ -f $READ1 && -f $READ2 ]]; then
            if [[ -f $LONG_READS ]]; then
                echo "Performing hybrid assembly with SPADES"
                spades.py -o ${OUTDIR}/01_assembly -1 ${READ1} -2 ${READ2} --${TECH} ${LONG_READS}
                exit_code=$?
                if [[ "$exit_code" -ne 0 ]]; then
                    echo "ERROR assembly with spades. (${exit_code})"
                    exit $exit_code
                fi 
            else
                echo "Performing short-read assembly with SPADES"
                spades.py -o ${OUTDIR}/01_assembly -1 ${READ1} -2 ${READ2}
                exit_code=$?
                if [[ "$exit_code" -ne 0 ]]; then
                    echo "ERROR assembly with spades. (${exit_code})"
                    exit $exit_code
                fi 
            fi
            # Rename SPADES' contigs output file name to the one expected by the pipeline.
            mv ${OUTDIR}/01_assembly/contigs.fasta ${contigs_fasta}
        else
            echo "Please check your input files."
            exit 1
        fi
        ;;
    "flye")
        case $TECH in
            "nanopore")
                FLYE_PARAMS="--nano-raw ${LONG_READS}"
                ;;
            "pacbio")
                FLYE_PARAMS="--pacbio-raw ${LONG_READS}"
                ;;
            *)
                echo "Unknown long read Technology -tech ${TECH}. Please use \"nanopore\" / \"pacbio\""
                ;;
        esac
        if [[ -n "$GENOME_SIZE" ]]; then
            FLYE_PARAMS="--genome-size ${GENOME_SIZE} ${FLYE_PARAMS}"
        fi
        flye --out-dir ${OUTDIR}/01_assembly --threads ${THREADS} ${FLYE_PARAMS}
        exit_code=$?
        if [[ "$exit_code" -ne 0 ]]; then
            echo "ERROR assembly with Flye. (${exit_code})"
            exit $exit_code
        fi 
        # Rename Flye's contigs output file name to the one expected by the pipeline.
        mv ${OUTDIR}/01_assembly/assembly.fasta ${contigs_fasta}
        ;;
    "skip")
        if [[ ! -f ${contigs_fasta} ]]; then
            echo "Cannot skip assembly as ${contigs_fasta} does not exists."
            exit 1
        fi
        ;;
    *)
        echo "Unknown assembler ${asm}. Please use \"flye\" / \"canu\" / \"wtdbg2\" / \"spades\""
        usage
        ;;
esac

# Racon Polishing
mkdir -p ${OUTDIR}/02_polish/01_racon
for ((i=0; i<=RACON_ITER; i++)); do
    if [[ "$i" -gt 0 ]]; then
        overlap_file="${OUTDIR}/02_polish/01_racon/overlap_iter_${i}.paf"
        minimap2 \
            -x ${MINIMAP2_PRESET} \
            -o ${overlap_file} \
            ${contigs_fasta} ${LONG_READS}
        exit_code=0 #$?
        if [[ "$exit_code" -ne 0 ]]; then
            break
        fi
        polished_sequence="${OUTDIR}/02_polish/01_racon/racon_iter_${i}.fa"
        racon \
            ${LONG_READS} \
            ${overlap_file} \
            ${contigs_fasta} > ${polished_sequence}
        exit_code=0 #$?
        if [[ "$exit_code" -ne 0 ]]; then
            echo "ERROR polishing with racon."
            break
        else
            contigs_fasta="${polished_sequence}"
        fi
    fi
done

if [[ $PILON_ITER -gt 0 || -f $READ1 || -f $READ2 ]]; then
    mkdir -p ${OUTDIR}/02_polish/02_pilon
    for ((i=0; i<=PILON_ITER; i++)); do
        if [[ $i -gt 0 ]]; then
            alignment_file="${OUTDIR}/02_polish/02_pilon/iter_${i}/alignment.bam"
            mkdir ${OUTDIR}/02_polish/02_pilon/iter_${i}

            bwa-mem2 index -p ${OUTDIR}/02_polish/02_pilon/iter_${i}/bwa_idx ${contigs_fasta}
            bwa-mem2 mem -t $HALF_THREADS ${OUTDIR}/02_polish/02_pilon/iter_${i}/bwa_idx $READ1 $READ2 | samtools sort -@ $((HALF_THREADS - 1)) > ${alignment_file}
            samtools index ${alignment_file}
            exit_code=$?
            if [[ "$exit_code" -ne 0 ]]; then
                echo "ERROR alignment with bwa-mem2."
                break
            fi
            pilon \
                -Xmx24G \
                --genome ${contigs_fasta} \
                --output pilon_polished \
                --outdir ${OUTDIR}/02_polish/02_pilon/iter_${i} \
                --bam ${alignment_file} \
                --changes --iupac 2>&1 | tee ${OUTDIR}/02_polish/02_pilon/iter_${i}/pilon.log
            exit_code=$?
            if [[ "$exit_code" -ne 0 ]]; then
                echo "ERROR polishing with pilon."
                break
            else
                polished_sequence="${OUTDIR}/02_polish/02_pilon/iter_${i}/pilon_polished.fasta"
                contigs_fasta="${polished_sequence}"
            fi
        fi
    done
fi
