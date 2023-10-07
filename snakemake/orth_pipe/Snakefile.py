import os
from scripts import printMetaDataInfo

PIPELINE_VERSION = "v2.1.9"

OUTPUT_FOLDER = os.path.abspath(config["OUTPUT_FOLDER"])

OGVERSION = config["OGVERSION"]

PREPINPUT_FOLDER = OUTPUT_FOLDER+"/"+config["PREP_FOLDER"]
ORTHO_FOLDER = OUTPUT_FOLDER+"/"+config["ORTHO_FOLDER"]
ORTHOFINDER_RESULTS = ORTHO_FOLDER+"/OrthoFinder"
SPECIES_TREE = ORTHOFINDER_RESULTS+"/Species_Tree/SpeciesTree_rooted_node_labels.txt"
#SPECIES_TREE = "/usr_storage/jcf/2.zhugecai/00.diploid/14.gene_family/3.gene_family/2.ortho//orthofinder/SpeciesTree_rooted_node_labels.txt"
TREEBEST = ORTHO_FOLDER+"/"+config["TREEBEST_FOLDER"]
#IADHORE_RESULTS = OUTPUT_FOLDER+"/"+config["SYNTENY_FOLDER"]+"/Results"
ALL_SPECIES = [key for key, value in config["inputFiles"].items()]
if "preWGD" in config:
  PREWGD = config["preWGD"]
else:
  PREWGD = []
if "postWGD" in config:
  POSTWGD = config["postWGD"]
else:
  POSTWGD= []
  if (config["FINDDUP"] == True):
    config["FINDDUP"] = False
if (type(PREWGD) != list):
  PREWGD = [PREWGD]
if (type(POSTWGD) != list):
  POSTWGD = [POSTWGD]
if len(POSTWGD)==1 or len(PREWGD)==0:
  config["FINDDUP"] = False
CLANS = [*PREWGD,*POSTWGD]

MACSE_PROGRAM=config["MACSE_PROGRAM"]
TREEBEST_PROGRAM=config["TREEBEST_PROGRAM"]
ORTHOFINDER_PROGRAM=config["ORTHOFINDER_PROGRAM"]
MAFFT_PROGRAM=config["MAFFT_PROGRAM"]
ALIGNMENT=ORTHO_FOLDER+"/"+config["ALIGN_FOLDER"]

NBTHREADS_IADHORE=config["NBTHREADS_IADHORE"]
IADHORE_PROGRAM = config["I-ADHORE_PROGRAM"]
IADHORE_INPUT = OUTPUT_FOLDER+"/"+config["SYNTENY_FOLDER"]+"/Input"

CDS_TREES = ORTHO_FOLDER+"/Treebest/allTrees.RDS"
CLANS_TREES = [ORTHO_FOLDER+"/Treebest/trees_dup_clans",ORTHO_FOLDER+"/Treebest/trees_nondup_clans"] if (config["FINDDUP"]) else [ORTHO_FOLDER+"/Treebest/trees_clans"]
TREES = CLANS_TREES if (config["EXTRACT_CLANS"]) else CDS_TREES

if os.path.exists(OUTPUT_FOLDER+"/metaData.txt"):
  os.remove(OUTPUT_FOLDER+"/metaData.txt")

rule all:
    input:
        OUTPUT_FOLDER+"/metaData.txt",
        [PREPINPUT_FOLDER+"/cdsFastas/"+str(SPC)+".fna" for SPC in ALL_SPECIES],
        #[PREPINPUT_FOLDER+"/geneIDtbls/"+str(SPC)+"_IDtbl.tsv" for SPC in ALL_SPECIES],
        #[PREPINPUT_FOLDER+"/genePosTbls/"+str(SPC)+"_genePos.tsv" for SPC in ALL_SPECIES],
        [PREPINPUT_FOLDER+"/proteinFastas/"+str(SPC)+".faa" for SPC in ALL_SPECIES],
        TREES,
        ORTHO_FOLDER+"/OGtbl.tsv",
        #IADHORE_RESULTS+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/Ohnolog_Tbl.tsv",
        #IADHORE_RESULTS+"/"+config["SYNTENY_CLANS_FOLDER"]+"/Ohnolog_Tbl.tsv" 

   
rule printMetaInfo:
  params:
    config = config,
    verion = PIPELINE_VERSION
  output:
    OUTPUT_FOLDER+"/metaData.txt"
  script:
    "scripts/printMetaDataInfo.R"

rule prepare_input:
   params:
      spc = lambda w: '{}'.format(w.SPC),
      protein_fasta = lambda w: config["inputFiles"]['{}'.format(w.SPC)]['protein_fasta'],
      cds_fasta = lambda w: config["inputFiles"]['{}'.format(w.SPC)]['cds_fasta'],
      gff = lambda w: config["inputFiles"]['{}'.format(w.SPC)]['gff'],
      source = lambda w: config["inputFiles"]['{}'.format(w.SPC)]['source']
   output:
      PREPINPUT_FOLDER+"/cdsFastas/{SPC}.fna",
      #PREPINPUT_FOLDER+"/geneIDtbls/{SPC}_IDtbl.tsv",
      #PREPINPUT_FOLDER+"/genePosTbls/{SPC}_genePos.tsv",
      PREPINPUT_FOLDER+"/proteinFastas/{SPC}.faa"
   script:
      "scripts/prepInputFiles.R"
        
#rule run_orthofinder:
#  input:
#        [PREPINPUT_FOLDER+"/proteinFastas/"+str(SPC)+".faa" for SPC in ALL_SPECIES]
#  params:
#        nbthreads = config["NBTHREADS_ORTHOFINDER"],
#        tmpdir = config["TEMP_FOLDER"]
#  output:
#        orthofinder_og_table = ORTHOFINDER_RESULTS+"/Orthogroups/Orthogroups.tsv",
#        dup_events = ORTHOFINDER_RESULTS+"/Gene_Duplication_Events/Duplications.tsv",
#        resolved_trees = directory(ORTHOFINDER_RESULTS+"/Resolved_Gene_Trees"),
#        species_tree = SPECIES_TREE
#  threads:
#        config["NBTHREADS_ORTHOFINDER"]
#  shell:
#        """
#        cp /usr_storage/jcf/2.zhugecai/00.diploid/14.gene_family/3.gene_family/2.ortho//orthofinder/pep/OrthoFinder/Results_Mar28//Orthogroups/Orthogroups.tsv {output.orthofinder_og_table}
#        cp /usr_storage/jcf/2.zhugecai/00.diploid/14.gene_family/3.gene_family/2.ortho//orthofinder/pep/OrthoFinder/Results_Mar28//Gene_Duplication_Events/Duplications.tsv {output.dup_events}
#        cp /usr_storage/jcf/2.zhugecai/00.diploid/14.gene_family/3.gene_family/2.ortho//orthofinder/pep/OrthoFinder/Results_Mar28///Resolved_Gene_Trees  {output.resolved_trees} -r
#        """
        
rule add_OG_version:
  input:
        orthofinder_og_table = ORTHOFINDER_RESULTS+"/Orthogroups/Orthogroups.tsv"
  output:
        OG_table_with_version = ORTHOFINDER_RESULTS+"/Orthogroups/Orthogroups_version.tsv"
  params:
        version = OGVERSION
  script:
        "scripts/add_version.R"
        
checkpoint extract_cds_orthogroups:
  input:
        OG_table_with_version = ORTHOFINDER_RESULTS+"/Orthogroups/Orthogroups_version.tsv",
        species_tree = SPECIES_TREE,
        dup_events = ORTHOFINDER_RESULTS+"/Gene_Duplication_Events/Duplications.tsv",
        resolved_trees = ORTHOFINDER_RESULTS+"/Resolved_Gene_Trees",
        cds_fasta = [PREPINPUT_FOLDER+"/cdsFastas/"+str(SPC)+".fna" for SPC in ALL_SPECIES]
  params:
        splitBefore = config["SPLIT_OG_BEFORE_CDS_ALIGNMENT"],
        support = config["SUPPORT_CLANS_TREES"],
        version = OGVERSION,
        clans_cds = config["CDS_ALIGNMENT_CLANS_ONLY"],
        max_nbSeq_macse = config["MAX_NBSEQ_MACSE"],
        max_lenSeq_macse = config["MAX_LENSEQ_MACSE"],
        preWGD = PREWGD,
        postWGD = POSTWGD,
        inputFiles = config["inputFiles"]
  output:
        cds_orthogroup_macse = directory(ALIGNMENT+"/cds_orthogroups_macse"),
        cds_orthogroup_mafft = directory(ALIGNMENT+"/cds_orthogroups_mafft")
  params:
        max_nbSeq_macse = config["MAX_NBSEQ_MACSE"],
        max_lenSeq_macse = config["MAX_LENSEQ_MACSE"]
  script:
        "scripts/extract_cds_orthogroups.R"
        
rule cds_alignment_macse:
  input:
    cds_orthogroup = ALIGNMENT+"/cds_orthogroups_macse/{OG}.fasta"
  output:
    nt_aln = ALIGNMENT+"/nt_aln/{OG}.fasta",
    aa_aln = ALIGNMENT+"/aa_aln/{OG}.fasta",
    nt_aln_tree_best = ALIGNMENT+"/nt_aln_for_treebest/{OG}.fasta",
  log:
    ALIGNMENT+"/macse_log/{OG}.log",
  shell:
        """
        set +e
        echo "" >> {log}
        echo "$(date)  -- Running macse on {input}" >> {log}
        echo "PID: $$" >> {log}
        
        java -jar {MACSE_PROGRAM} -prog alignSequences -seq {input.cds_orthogroup} -out_NT {output.nt_aln} -out_AA {output.aa_aln} >> {log} 2>&1
        exitcode=$?
        # increase memory if it fails:
        if [ $exitcode -ne 0 ] || [ ! -f {output.nt_aln} ]
        then
          echo "Error: MACSE failed. trying again with more memory" >> {log}
          java -Xmx8000m -jar {MACSE_PROGRAM} -prog alignSequences -seq {input.cds_orthogroup} -out_NT {output.nt_aln} -out_AA {output.aa_aln} >> {log} 2>&1
        fi
        Rscript scripts/removeFrameShiftStopCodon.R {output.nt_aln} {output.aa_aln} {output.nt_aln_tree_best} >> {log} 2>&1
        """
        
rule cds_alignment_mafft:
  input:
        cds_orthogroup = ALIGNMENT+"/cds_orthogroups_mafft/{OG}.fasta"
  output:
        nt_aln = ALIGNMENT+"/nt_aln/{OG}.fasta",
        aa_aln = ALIGNMENT+"/aa_aln/{OG}.fasta",
        nt_aln_tree_best = ALIGNMENT+"/nt_aln_for_treebest/{OG}.fasta",
  shell:
        """
        java -jar {MACSE_PROGRAM} -prog translateNT2AA -seq {input.cds_orthogroup} -out_AA {ALIGNMENT}/aa_aln/{wildcards.OG}_tmp.fasta
        {MAFFT_PROGRAM} --quiet {ALIGNMENT}/aa_aln/{wildcards.OG}_tmp.fasta > {output.aa_aln}
        rm {ALIGNMENT}/aa_aln/{wildcards.OG}_tmp.fasta
        {TREEBEST_PROGRAM} backtrans {output.aa_aln} {input.cds_orthogroup} > {output.nt_aln_tree_best}
        cp {output.nt_aln_tree_best} {output.nt_aln}
        """
        
rule run_treebest:
  input:
        nt_aln_tree_best = ALIGNMENT+"/nt_aln_for_treebest/{OG}.fasta",
        species_tree = SPECIES_TREE
  output:
        tree_best = TREEBEST+"/trees/{OG}.nhx",
        distmat = TREEBEST+"/distmat/{OG}.dist"
  shell:
        """
        set +e
        {TREEBEST_PROGRAM} distmat ds {input.nt_aln_tree_best} > {output.distmat}
        {TREEBEST_PROGRAM} best -f {input.species_tree} -o {output.tree_best} {input.nt_aln_tree_best}
        exitcode=$?
        if [ $exitcode -ne 0 ] || [ ! -f {output.tree_best} ]
        then
            {TREEBEST_PROGRAM} best -F 1 -f {input.species_tree} -o {output.tree_best} {input.nt_aln_tree_best}
            exitcode=$?
            if [ $exitcode -ne 0 ] || [ ! -f {output.tree_best} ]
            then
               {TREEBEST_PROGRAM} best -Z 0.0001 -F 1 -f {input.species_tree} -o {output.tree_best} {input.nt_aln_tree_best}
            fi
        fi
        exit 0
        """
        
def aggregate_input(s):
    out_macse = checkpoints.extract_cds_orthogroups.get().output[0]
    out_mafft = checkpoints.extract_cds_orthogroups.get().output[1]
    return(expand(TREEBEST+"/trees/{OG}.nhx",OG=
        [*glob_wildcards(os.path.join(out_macse, "{OG}.fasta")).OG, *glob_wildcards(os.path.join(out_mafft, "{OG}.fasta")).OG]))
        
rule all_trees_in_one_file:
  input:
        aggregate_input
  output:
        nhx_tree = TREEBEST+"/alltrees_nhx.tsv",
        rds_tree = TREEBEST+"/allTrees.RDS"
  script:
        "scripts/all_trees_in_one_file.R"
 
rule split_OG_table:
  input:
        rds_tree = TREEBEST+"/allTrees.RDS",
        species_tree = SPECIES_TREE
  output:
        og_tbl_file_out = ORTHO_FOLDER+"/OGtbl.tsv"
  script:
        "scripts/makeSplitOGtable.R"

rule extract_clans_dup:
  input:
        rds_tree = TREEBEST+"/allTrees.RDS",
        species_tree = SPECIES_TREE,
        og_tbl = ORTHO_FOLDER+"/OGtbl.tsv"
  params:
        supportClan = config["SUPPORT_CLANS_TREES"],
        supportDup = config["SUPPORT_DUP_TREES"],
        min = config["MIN_CLANS_SIZE"],
        preWGD = PREWGD,
        postWGD = POSTWGD
  output:
        dup_clans_tree_nwk = directory(TREEBEST+"/trees_dup_clans"),
        nondup_clans_tree_nwk = directory(TREEBEST+"/trees_nondup_clans"),
        dup_clans_cds = directory(ALIGNMENT+"/nt_aln_dup_clans"),
        nondup_clans_cds = directory(ALIGNMENT+"/nt_aln_nondup_clans")
  script:
        "scripts/extract_clans_dup.R"

rule extract_clans:
  input:
        rds_tree = TREEBEST+"/allTrees.RDS",
        species_tree = SPECIES_TREE,
        og_tbl = ORTHO_FOLDER+"/OGtbl.tsv"
  params:
        supportClan = config["SUPPORT_CLANS_TREES"],
        min = config["MIN_CLANS_SIZE"],
        preWGD = PREWGD,
        postWGD = POSTWGD
  output:
        clans_tree_nwk = directory(TREEBEST+"/trees_clans"),
        clans_cds = directory(ALIGNMENT+"/nt_aln_clans")
  script:
        "scripts/extract_clans.R"

        
#rule iadhore_input_whole_OG:
#  input:
#        orthofinder_og_table = ORTHOFINDER_RESULTS+"/Orthogroups/Orthogroups_version.tsv",
#        genePos = [PREPINPUT_FOLDER+"/genePosTbls/"+str(SPC)+"_genePos.tsv" for SPC in CLANS]
#  output:
#        geneListDir = [directory(IADHORE_INPUT+"/geneLists/"+str(SPC)) for SPC in CLANS],
#        gene_family = IADHORE_INPUT+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/gene_family.txt",
#        setting = IADHORE_INPUT+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/settings.ini"
#  script:
#        "scripts/prepIADHORE_whole_OG.R"

#rule run_iadhore_whole_OG:
#  input:
#        geneListDir = [IADHORE_INPUT+"/geneLists/"+str(SPC) for SPC in CLANS],
#        gene_family = IADHORE_INPUT+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/gene_family.txt",
#        setting = IADHORE_INPUT+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/settings.ini"
#  output:
#        anchorpoints = IADHORE_RESULTS+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/anchorpoints.txt",
#        multiplicons = IADHORE_RESULTS+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/multiplicons.txt",
#        genes = IADHORE_RESULTS+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/genes.txt"
#  params:
#        MPI = config["IADHORE_MPI"]
#  threads:
#        NBTHREADS_IADHORE
#  shell:
#        """
#        scripts/run_iadhore.sh {threads} {IADHORE_PROGRAM} {input.setting} {params.MPI}
#        """        

#rule ohnolog_table:
#  input:
#      anchorpoints = IADHORE_RESULTS+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/anchorpoints.txt",
#      multiplicons = IADHORE_RESULTS+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/multiplicons.txt",
#      genes = IADHORE_RESULTS+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/genes.txt",
#      og_tbl = ORTHO_FOLDER+"/OGtbl.tsv",
#      rds_tree = TREEBEST+"/allTrees.RDS"
#  params:
#      level2 = config["LEVEL2_ONLY"]
#  output:
#      ohnolog_tbl = IADHORE_RESULTS+"/"+config["SYNTENY_WHOLE_OG_FOLDER"]+"/Ohnolog_Tbl.tsv" 
#  script:
#      "scripts/generate_ohnolog_tbl.R"
      
      
#rule iadhore_input_clans:
#  input:
#        species_tree = SPECIES_TREE,
#        og_tbl = ORTHO_FOLDER+"/OGtbl.tsv",
#        geneListDir = [IADHORE_INPUT+"/geneLists/"+str(SPC) for SPC in CLANS],
#        genePos = [PREPINPUT_FOLDER+"/genePosTbls/"+str(SPC)+"_genePos.tsv" for SPC in CLANS]
#  output:
#        gene_family = IADHORE_INPUT+"/"+config["SYNTENY_CLANS_FOLDER"]+"/gene_family.txt",
#        setting = IADHORE_INPUT+"/"+config["SYNTENY_CLANS_FOLDER"]+"/settings.ini"
#  script:
#        "scripts/prepIADHORE_clans.R"


#rule run_iadhore_clans:
#  input:
#        geneListDir = [IADHORE_INPUT+"/geneLists/"+str(SPC) for SPC in CLANS],
#        gene_family = IADHORE_INPUT+"/"+config["SYNTENY_CLANS_FOLDER"]+"/gene_family.txt",
#        setting = IADHORE_INPUT+"/"+config["SYNTENY_CLANS_FOLDER"]+"/settings.ini"
#  output:
#        anchorpoints = IADHORE_RESULTS+"/"+config["SYNTENY_CLANS_FOLDER"]+"/anchorpoints.txt",
#        multiplicons = IADHORE_RESULTS+"/"+config["SYNTENY_CLANS_FOLDER"]+"/multiplicons.txt",
#        genes = IADHORE_RESULTS+"/"+config["SYNTENY_CLANS_FOLDER"]+"/genes.txt"
#  threads:
#        NBTHREADS_IADHORE
#  shell:
#        """
#        scripts/run_iadhore.sh {threads} {IADHORE_PROGRAM} {input.setting}
#        """
      
#rule ohnolog_table_clans:
#  input:
#      anchorpoints = IADHORE_RESULTS+"/"+config["SYNTENY_CLANS_FOLDER"]+"/anchorpoints.txt",
#      multiplicons = IADHORE_RESULTS+"/"+config["SYNTENY_CLANS_FOLDER"]+"/multiplicons.txt",
#      genes = IADHORE_RESULTS+"/"+config["SYNTENY_CLANS_FOLDER"]+"/genes.txt",
#      og_tbl = ORTHO_FOLDER+"/OGtbl.tsv",
#      rds_tree = TREEBEST+"/allTrees.RDS"
#  params:
#      level2 = config["LEVEL2_ONLY"]
#  output:
#      ohnolog_tbl = IADHORE_RESULTS+"/"+config["SYNTENY_CLANS_FOLDER"]+"/Ohnolog_Tbl.tsv"
#  script:
#      "scripts/generate_ohnolog_tbl.R"
           

        
        
