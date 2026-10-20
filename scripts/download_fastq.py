
#!/usr/bin/env python3

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


PROJECT = Path(__file__).resolve().parent.parent
SAMPLESHEET = PROJECT / "data" / "metadata" / "shotgun_samples.tsv"
RAW_DIR = PROJECT / "data" / "raw"

# External disk: temporary files and SRA cache
EXTERNAL_BASE = Path(
    "/Volumes/Paul/LARISA/github_projects"
)
DEFAULT_TEMP_DIR = (
    Path("/Volumes/Paul/LARISA/github_projects/shotgun_metagenomics_work")
    / "sra_download_tmp"
)
DEFAULT_CACHE_DIR = EXTERNAL_BASE / "sra_cache"


def run_command(command):
    print("+", " ".join(map(str, command)), flush=True)
    subprocess.run(command, check=True)


def valid_fastq_pair(read1, read2):
    """Check gzip integrity, FASTQ structure and matching read counts."""
    if not read1.is_file() or not read2.is_file():
        return False

    if read1.stat().st_size == 0 or read2.stat().st_size == 0:
        return False

    counts = []

    try:
        for path in (read1, read2):
            with gzip.open(path, "rt") as handle:
                lines = sum(1 for _ in handle)

            if lines == 0 or lines % 4 != 0:
                return False

            counts.append(lines // 4)

        return counts[0] == counts[1]

    except (OSError, EOFError, gzip.BadGzipFile):
        return False


def download_sample(run, sample, raw_dir, cache_dir, temp_base, threads):
    read1 = raw_dir / f"{sample}_1.fastq.gz"
    read2 = raw_dir / f"{sample}_2.fastq.gz"

    if read1.exists() or read2.exists():
        if valid_fastq_pair(read1, read2):
            print(f"[SKIP] {sample}: FASTQ completos y validados.")
            return True

        print(
            f"[ERROR] {sample}: FASTQ existentes incompletos o inválidos. "
            "No se borrarán automáticamente.",
            file=sys.stderr,
        )
        return False

    print(f"\n[START] {sample} ({run})", flush=True)

    cache_dir.mkdir(parents=True, exist_ok=True)
    temp_base.mkdir(parents=True, exist_ok=True)

    try:
        run_command([
            "prefetch", run,
            "--output-directory", str(cache_dir),
        ])

        sra_file = cache_dir / run / f"{run}.sra"

        if not sra_file.is_file():
            raise RuntimeError(f"No se encontró el archivo SRA: {sra_file}")

        with tempfile.TemporaryDirectory(
            prefix=f"{sample}_", dir=temp_base
        ) as temp_name:
            temp_dir = Path(temp_name)

            run_command([
                "fasterq-dump", str(sra_file),
                "--split-files",
                "--threads", str(threads),
                "--temp", str(temp_dir),
                "--outdir", str(temp_dir),
            ])

            temp_r1 = temp_dir / f"{run}_1.fastq"
            temp_r2 = temp_dir / f"{run}_2.fastq"

            if not temp_r1.is_file() or not temp_r2.is_file():
                raise RuntimeError(
                    "No se generaron ambos FASTQ pareados. "
                    "Comprueba si la accesión es paired-end."
                )

            temp_gz1 = temp_dir / f"{sample}_1.fastq.gz"
            temp_gz2 = temp_dir / f"{sample}_2.fastq.gz"

            for source, destination in (
                (temp_r1, temp_gz1),
                (temp_r2, temp_gz2),
            ):
                with open(source, "rb") as src:
                    with gzip.open(destination, "wb") as dst:
                        shutil.copyfileobj(src, dst)

            if not valid_fastq_pair(temp_gz1, temp_gz2):
                raise RuntimeError(
                    "Los FASTQ no pasaron la validación de integridad "
                    "y número de lecturas."
                )

            # Publish final files only after validation.
            shutil.move(str(temp_gz1), str(read1))

            try:
                shutil.move(str(temp_gz2), str(read2))
            except Exception:
                read1.unlink(missing_ok=True)
                raise

        print(f"[OK] {sample}: descarga validada.")
        return True

    except (subprocess.CalledProcessError, OSError, RuntimeError) as exc:
        print(f"[ERROR] {sample}: {exc}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Descarga secuencial de FASTQ desde SRA."
    )
    parser.add_argument(
        "--tmp-dir",
        type=Path,
        default=DEFAULT_TEMP_DIR,
        help="Directorio temporal en el disco externo.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help="Caché SRA en el disco externo.",
    )
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()

    if not SAMPLESHEET.is_file():
        sys.exit(f"No existe la tabla: {SAMPLESHEET}")

    if not Path("/Volumes/Paul").is_dir():
        sys.exit(
            "No se encuentra /Volumes/Paul. "
            "Conecta el disco externo antes de descargar."
        )

    for executable in ("prefetch", "fasterq-dump"):
        if shutil.which(executable) is None:
            sys.exit(
                f"No se encontró {executable}. "
                "Activa el entorno que contiene SRA Toolkit."
            )

    RAW_DIR.mkdir(parents=True, exist_ok=True)

    with SAMPLESHEET.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "participant", "group", "timepoint", "run"}

        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            sys.exit(
                "La tabla debe contener: "
                "sample, participant, group, timepoint, run."
            )

        samples = list(reader)

    if not samples:
        sys.exit("La tabla no contiene muestras.")

    sample_names = [row["sample"].strip() for row in samples]
    if len(sample_names) != len(set(sample_names)):
        sys.exit("Hay nombres de muestra duplicados en la tabla.")

    failures = []

    for row in samples:
        sample = row["sample"].strip()
        run = row["run"].strip()

        if not sample or not run:
            print(f"[ERROR] Fila incompleta: {row}", file=sys.stderr)
            failures.append(sample or str(row))
            continue

        success = download_sample(
            run=run,
            sample=sample,
            raw_dir=RAW_DIR,
            cache_dir=args.cache_dir,
            temp_base=args.tmp_dir,
            threads=args.threads,
        )

        if not success:
            failures.append(sample)

    print("\n=== RESUMEN ===")
    print(f"Muestras en la tabla: {len(samples)}")
    print(f"Fallos: {len(failures)}")

    if failures:
        print("Revisar:", ", ".join(failures))
        sys.exit(1)

    print("Todas las muestras están descargadas y validadas.")


if __name__ == "__main__":
    main()

