import sys
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt
import shutil

from ancestry_to_23andme import process
from plotter import plotter

#Verbose / -v flag to control the amount of info printed in std output.
VERBOSE = False


def find_program(display_name, possible_names):
    """
    Setting up a way to find PLINK1.9 and PLINK2 since the name
    of installed plink tools are not the same on all machines.
    This is a helper function used later to check dependencies
    to ensure that both PLINK versions are installed and
    accessible.
    """
    for name in possible_names:
        path = shutil.which(name)

        if path is not None:
            return path

    #Names to display if either PLINK version is not found
    names = ", ".join(possible_names)
    raise FileNotFoundError(
        f"{display_name} was not found.\n"
        f"Expected one of these commands: {names}\n"
        "Install it and ensure it is available on your PATH."
    )


def check_dependencies():
    """
    Uses the find_program helper to find appropriate versions of PLINK
    and assign each of them to one name that can be accessed by other
    program/function calls later on.
    """
    plink19 = find_program(
        "PLINK 1.9",
        ["plink1.9", "plink"]
    )

    plink2 = find_program(
        "PLINK 2",
        ["plink2"]
    )

    return plink19, plink2






def parse_args():
    args = sys.argv[1:]

    verbose = "-v" in args or "--verbose" in args

    args = [arg for arg in args if arg not in ('-v', '--verbose')]

    if "-o" not in args:
        print("Usage:")
        print("\tpython3 run_relatedness.py file1.txt file2.txt [file3.txt ...] -o output_folder")
        print("OR")
        print("\tpython3 run_relatedness.py folder_with_dna_files -o output_folder")
        print()
        print(f"Note: Optionally use '-v' or '--verbose' to display detailed PLINK output.")
        sys.exit(1)


    o_index = args.index("-o")
    input_args = args[:o_index]

    if o_index + 1 >= len(args):
        print("Error: -o was given but no output folder name was provided.")
        sys.exit(1)

    output_name = args[o_index + 1]

    if not input_args:
        print("Error: Must provide DNA files (as txt) or a folder before -o.")
        sys.exit(1)

    return input_args, output_name, verbose


def run(command, description, quiet_failure=False):
    """
    An easy way for me to run PLINK CLI commands.
    quiet_failure is intended for the PLINK merge where one error is often times expected in the case of poly-allelic SNPs.
    """
    print(f"\n{description}...")

    if VERBOSE:
        print("Running:")
        print(" ".join(map(str, command)))
        subprocess.run(command, check=True)

    else:
        try:
            subprocess.run(command, check=True, text=True, stdout=subprocess.DEVNULL, stderr=(subprocess.DEVNULL if quiet_failure else None))

        except subprocess.CalledProcessError as error:
            if not quiet_failure:
                print(f"\nError: <{description}> failed.")

            raise 


#Returns a list of files corresponding to each DNA file intended for pairwise analysis
def get_DNA(args):
    paths = [Path(i) for i in args]

    #IF we're working with ONE folder:
    if len(paths) == 1 and paths[0].exists() and paths[0].is_dir():
        folder = paths[0]

        dna_files = [i for i in folder.iterdir() if i.is_file() and i.suffix.lower() in {".txt", ".tsv"}]

        if len(dna_files) < 2:
            print("Need at least 2 DNA files for pairwise analysis. Please check your folder and retry.")
            sys.exit(1)

        print()
        print(f"Found {len(dna_files)} possible DNA files in folder '{folder}':")
        for f in dna_files:
            print("\t-", f.name)

        return dna_files
    

    #IF instead we're working with individual files:
    dna_files = []

    for path in paths:
        if not path.exists():
            print(f"Could not find input: {path}")
            sys.exit(1)

        if not path.is_file():
            print(f"Input is not a file: {path}")
            print("Provide either one folder, or a series of DNA files.")
            sys.exit(1)

        if path.suffix.lower() not in {".txt", ".tsv"}:
            print(f"Unsupported file type: {path}")
            print("Allowed file types: .txt and .tsv")
            sys.exit(1)

        dna_files.append(path)

    if len(dna_files) < 2:
        print("Need at least 2 DNA files for pairwise analysis.")
        sys.exit(1)

    print()
    print(f"Using {len(dna_files)} DNA files:")
    for f in dna_files:
        print("\t-", f.name)

    return dna_files


#Stores the name of each file for stable PLINK output. Ensures name formatting is usable.
def conserve_name(path):
    name = path.stem
    name = name.replace(" ", "_")
    name = "".join(c for c in name if c.isalnum() or c in {"_", "-"})
    return name


#Returns a list of 2-tuples where tup[0] = sample_name and tup[1] = File
def convert_all_DNA(dna_files, outdir):
    """
    Convert each supplied DNA file into a PLINK-compatible text file.
    Returns a list of tuples:
        (sample_name, converted_file)
    """

    converted_root = outdir / "converted_rawDNA"
    converted_root.mkdir(parents=True, exist_ok=True)

    converted_files = []

    for dna_file in dna_files:
        sample_name = conserve_name(dna_file)

        sample_converted_dir = converted_root / sample_name
        sample_converted_dir.mkdir(parents=True, exist_ok=True)

        # process() adds ".txt" itself
        converted_prefix = sample_converted_dir / sample_name

        print()
        print(f"Converting: {dna_file.name}")

        converted_txt = process(
            dna_file,
            converted_prefix
        )

        # process() currently returns a string, so convert it to Path
        converted_txt = Path(converted_txt)
        if not converted_txt.exists() or converted_txt.stat().st_size == 0:
            print()
            print(f"Error: No valid SNP records were recognised in '{dna_file.name}'.")
            print("The file may use an unsupported raw-DNA format, or not contain any usable information.")
            print("Check README for information on supported file types.")
            sys.exit(1)

        print("Converted file saved to:", converted_txt)

        converted_files.append(
            (sample_name, converted_txt)
        )

    return converted_files


#Uses tup[1] to create plink (bed, bam, fam) files for each sample.
def make_plink_files(converted_files, outdir):
    plink_root = outdir / "plink_samples"
    plink_root.mkdir(parents=True, exist_ok=True)

    plink_prefixes = []

    for sample_name, converted_DNA in converted_files:

        # Each sample gets its own PLINK folder
        sample_plink_dir = plink_root / sample_name
        sample_plink_dir.mkdir(parents=True, exist_ok=True)

        plink_prefix = sample_plink_dir / sample_name

        run([
            PLINK19,
            "--23file", str(converted_DNA), sample_name, sample_name,
            "--snps-only", "just-acgt",
            "--autosome",
            "--make-bed",
            "--out", str(plink_prefix)
        ], description=f"Creating PLINK files for {sample_name}")

        plink_prefixes.append(plink_prefix)

    return plink_prefixes
    # A list of the three PLINK files stored for each sample/person


#Merge datasets. Uses missnp if merged SNP arrays are biallelic. Could use PLINK2.0 machinery but this is easier to follow and I understand it better + compatibility.
def merge(plink_files, outdir):

    merged_dir = outdir/"merged_files"
    merged_dir.mkdir(parents=True, exist_ok=True)

    base = plink_files[0]
    rest = plink_files[1:]

    merge_list = merged_dir / "merge_list.txt"

    with open(merge_list, "w") as f:
        for prefix in rest:
            f.write(f"{prefix}.bed {prefix}.bim {prefix}.fam\n")

    merged_prefix = merged_dir / "merged_all"


    try:
        run([
            PLINK19,
            "--bfile", str(base),
            "--merge-list", str(merge_list),
            "--make-bed",
            "--out", str(merged_prefix)
        ], description=f"Merging all genotype datasets", quiet_failure=True)

    except subprocess.CalledProcessError:
        print()
        print("Initial Merge Failed: Likely due to poly-allelic variants")
        print("Attempting Fallback: Exclude SNPs listed in the newly created .missnp file.")
        print()

        missnp_file = Path(str(merged_prefix) + "-merge.missnp")

        if not missnp_file.exists():
            print("Could not find:", missnp_file)
            print("The merge failed some other currently unresolvable reason. Contact me via my Github info if you desperately need a fix and I will try my best.")
            sys.exit(1)

        cleaned_root = outdir / "plink_samples_cleaned"
        cleaned_root.mkdir(parents=True, exist_ok=True)

        cleaned_prefixes = []

        for prefix in plink_files:
            sample_name = prefix.name

            sample_clean_dir = cleaned_root / sample_name
            sample_clean_dir.mkdir(parents=True, exist_ok=True)

            cleaned_prefix = sample_clean_dir / sample_name


            #Reconvert to PLINK binary formats using missnp
            run([
            PLINK19,
            "--bfile", str(prefix),
            "--exclude", str(missnp_file),
            "--make-bed",
            "--out", str(cleaned_prefix)
            ], description=f"Removing incompatible (poly-allelic) variants from {sample_name}")
            cleaned_prefixes.append(cleaned_prefix)

        base_clean = cleaned_prefixes[0]
        rest_clean = cleaned_prefixes[1:]

        cleaned_merge_list = merged_dir / "merge_list_cleaned.txt"

        with open(cleaned_merge_list, "w") as f:
            for prefix in rest_clean:
                f.write(f"{prefix}.bed {prefix}.bim {prefix}.fam\n")


        #Redo merge with triallelic + snps removed
        try:
            run([
            PLINK19,
            "--bfile", str(base_clean),
            "--merge-list", str(cleaned_merge_list),
            "--merge-equal-pos",
            "--make-bed",
            "--out", str(merged_prefix)
            ], description=f"Performing backup merge with poly-allelic variants removed\n")

        except subprocess.CalledProcessError:
            print("Second Merge Failed.")
            print("Post-poly-allele removal merge failed. The files are still incompatible for some reason (unsure yet - will fix later if this exception is raised)")
            sys.exit(1)

    return merged_prefix


def remove_duplicate_snps(merged_prefix, outdir):
    merged_no_dups = outdir / "merged_files" / "merged_all_no_dups"

    run([
    PLINK2,
    "--bfile", str(merged_prefix),
    "--set-all-var-ids", "@:#",
    "--rm-dup", "exclude-mismatch", "list",
    "--make-bed",
    "--out", str(merged_no_dups)
], description="Removing duplicate variants")

    return merged_no_dups


def run_king(merged_no_dups, outdir):
    relatedness_dir = outdir / "relatedness"
    relatedness_dir.mkdir(parents=True, exist_ok=True)

    relatedness_prefix = relatedness_dir / "relatedness_all"

    run([
        PLINK2,
        "--bfile", str(merged_no_dups),
        "--make-king-table",
        "--out", str(relatedness_prefix)
    ], description="Calculating pairwise relatedness with KING")

    kin_file = Path(str(relatedness_prefix) + ".kin0")

    print()
    print("KING results saved to:", kin_file)

    if not kin_file.exists():
        raise FileNotFoundError(f"KING did not produce the expected results file: {kin_file}")

    print()
    print("Raw KING relatedness results:")
    print(kin_file.read_text())

    return kin_file


def main():
    global VERBOSE
    global PLINK19
    global PLINK2

    input_args, output_name, VERBOSE = parse_args()
    PLINK19, PLINK2 = check_dependencies()
    
    dna_files = get_DNA(input_args)

    outdir = Path(output_name)
    outdir.mkdir(parents=True, exist_ok=True)

    converted_files = convert_all_DNA(dna_files, outdir)
    plink_files = make_plink_files(converted_files, outdir)

    merged_prefix = merge(plink_files, outdir)
    merged_no_dups = remove_duplicate_snps(merged_prefix, outdir)

    kin_file = run_king(merged_no_dups, outdir)

    user_df, fig = plotter(kin_file)

    csv_file = kin_file.parent / "pairwise_relatedness.csv"
    plot_file = kin_file.parent / "pairwise_relatedness.png"

    user_df.to_csv(csv_file, index=False)

    fig.savefig(plot_file, dpi=300, bbox_inches="tight")

    print()
    print("Readable relationship table saved to:", csv_file)
    print("Relatedness plot saved to:", plot_file)

    print()
    print("Done.")
    print("All final results saved in:", outdir)

    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    try:
        main()

    except subprocess.CalledProcessError:
        if VERBOSE:
            raise

        else:
            print()
            print("The analysis could not be completed.")
            print("Run again with --verbose to see detailed diagnostic output.")
            sys.exit(1)

    except FileNotFoundError as error:
        print()
        print(f"Error: {error}")
        sys.exit(1)

