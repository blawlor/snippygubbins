import argparse
from subprocess import PIPE
from subprocess import call
import sys
from os import listdir, removedirs, makedirs
from os.path import isfile, join, splitext, exists, basename
import shutil

def isIllumina(f):
    components = f.replace(".fastq.gz","").split("_")
    no_of_components = len(components)
    return no_of_components >= 5 and (components[no_of_components-2] == "R1" or components[no_of_components-2] == "R2")

def isOtherFormat(f):
    components = f.replace(".fastq.gz","").split("_")
    return len(components) == 2 and (components[1] == "1" or components[1] == "2")

def isValidFileFormat(f):
    if (isIllumina(f)):
        return True
    if (isOtherFormat(f)):
        return True
    return False

def isForwards(filename):
    components = filename.replace(".fastq.gz","").split("_")
    no_of_components = len(components)
    if (isIllumina(filename)):
        print(filename + " is illumina format")
        direction = components[no_of_components-2]
        print("Direction is " + direction)
        return direction == "R1"
    print(filename + " is other format")
    if (isOtherFormat(filename)):
        direction = components[1]
        print("Direction is " + direction)
        return direction == "1"
    return False


def isReverse(filename):
    components = filename.replace(".fastq.gz","").split("_")
    no_of_components = len(components)
    if (isIllumina(filename)):
        print(filename + " is illumina format")
        direction = components[no_of_components-2]
        print("Direction is " + direction)
        return direction == "R2"
    if (isOtherFormat(filename)):
        print(filename + " is other format")
        direction = components[1]
        print("Direction is " + direction)
        return direction == "2"
    return False

def extractSampleName(fastq_file):
    if (isIllumina(fastq_file)):
        components = fastq_file.replace(".fastq.gz", "").split("_")
        no_of_components = len(components)
        no_of_name_components = no_of_components - 4;
        name_components = components[0:no_of_name_components]
        return "_".join(name_components)
    else:
        return f.split("_")[0]

def createSamplesDict(fastq_files):
    #Create a dictionary of sample names and their files
    sample_prefixes = [extractSampleName(f) for f in fastq_files if isValidFileFormat(f)]
    result = {}
    for sample in sample_prefixes:
        print("Finding files for sample " + sample)
        r1_file = next(f for f in fastq_files if (f.startswith(sample) and isForwards(f)))
        print(str(r1_file))
        r2_file = next(f for f in fastq_files if (f.startswith(sample) and isReverse(f)))
        print(str(r2_file))
        result[sample] = [r1_file, r2_file]
    return result

def create_dir(directory):
    if not exists(directory):
        makedirs(directory)

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

def execute(indir, outdir, cpus, do_snippy, do_combinevcf, do_gubbins):
    print("Using " + str(cpus) + " cpus")
    #Call Snippy on each pair of fastq.gz paired end files in the input directory
    fastq_files = [f for f in listdir(indir) if (isfile(join(indir, f)) and f.lower().endswith(".fastq.gz"))]
    #Create a dictionary of sample names and paired end files.
    sample_dict = createSamplesDict(fastq_files)
    if do_snippy:
        print("About to perform Snippy alignment for all samples")
        for sample in sample_dict.keys():
            if not exists(outdir + "/" + sample):
                sample_outdir = outdir+"/"+sample+"_snippy"
                snippy_command = ["./bin/snippy", "--cpus", str(cpus), "--outdir", sample_outdir, "--ref", indir + "/ref.fasta",
                      "--R1", indir + "/" + sample_dict[sample][0],
                      "--R2", indir + "/" + sample_dict[sample][1]]
                call(snippy_command, stdout=sys.stdout, stderr=sys.stdout)
                print("Completed Snippy alignment of " + sample)
                #Copy the required output files to directories named per sample
                create_dir(outdir + "/" + sample)
                create_dir(outdir + "/" + sample+"/reference")
                shutil.move(sample_outdir+"/snps.vcf", outdir+"/"+sample+"/"+sample+".vcf")
                shutil.move(sample_outdir+"/snps.tab", outdir+"/"+sample+"/"+sample+".tab")
                shutil.move(sample_outdir+"/snps.aligned.fa", outdir+"/"+sample+"/"+sample+".aligned.fa")
                shutil.move(sample_outdir+"/reference/ref.fa", outdir+"/"+sample+"/reference/ref.fa")
                #Remove the large temp snippy files
                shutil.rmtree(outdir+"/"+sample+"_snippy")

    if do_combinevcf:
        # Combine the vcf files into a single multiVCF
        print("About to call snippy-core to create multivcf")
        snippy_core_command = ["./bin/snippy-core", "--prefix", "combined"]
        for sample in sample_dict.keys():
            snippy_core_command.append(outdir + "/" + sample)
        snippy_core_command_string = " ".join(str(x) for x in snippy_core_command)
        with open(outdir + "/snippycorecommand.txt", "w") as snippy_core_command_file:
            snippy_core_command_file.write(snippy_core_command_string)
        call(snippy_core_command, stdout=sys.stdout, stderr=sys.stdout)
        combined_files = [f for f in listdir(".") if (isfile(join(".", f)) and basename(f).lower().startswith("combined"))]
        for combined_file in combined_files:
            print("Copying " + combined_file + " to output dir " + outdir)
            shutil.copy(combined_file, outdir + "/" + combined_file)

        print("Completed multivcf generation for all samples")


    if do_gubbins:
        print("About to call gubbins")
        gubbins_command = ["run_gubbins",
                           "--tree_builder", "fasttree",
                             "--prefix", "gubbins",
                           outdir + "/combined.full.aln"
                           ]
        create_dir(outdir+"/tree")
        call(gubbins_command, stdout=sys.stdout, stderr=sys.stdout)
        gubbins_files = [f for f in listdir(".") if (isfile(join(".", f)) and basename(f).lower().startswith("gubbins."))]
        for gubbins_file in gubbins_files:
            shutil.move(gubbins_file, outdir + "/tree/" + gubbins_file)
        print("Created tree using Gubbins")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", type=str, help="Directory containing input files")
    parser.add_argument("--outdir", type=str, help="Directory (existing) in which output files will be placed")
    parser.add_argument("--cpus", type=int, help="Number of CPUs to be used by Snippy")
    parser.add_argument('--snippy', type=str, help="yes/true if you want snippy to be called")
    parser.set_defaults(snippy='yes')
    parser.add_argument('--combinevcf', type=str, help="yes/true if you want snippy core to be called")
    parser.set_defaults(combinevcf='yes')
    parser.add_argument('--gubbins', type=str, help="yes/true if you want to invoke gubbins")
    parser.set_defaults(gubbins='yes')
    args = parser.parse_args()
    execute(args.indir, args.outdir, args.cpus, str2bool(args.snippy), str2bool(args.combinevcf), str2bool(args.gubbins))
