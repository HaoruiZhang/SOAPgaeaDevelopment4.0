/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.mapreduce.haplotypecaller;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.bgi.flexlab.gaea.tools.genotyer.VariantCallingEngine;
import org.bgi.flexlab.gaea.tools.haplotypecaller.HaplotypeCallerEngine;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.GenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.SNPGenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.genotyer.genotypecaller.AFCalcFactory;
import org.bgi.flexlab.gaea.util.GaeaVariantContextUtils;
import org.bgi.flexlab.gaea.util.pairhmm.PairHMM;
import org.seqdoop.hadoop_bam.SAMFormat;
import org.seqdoop.hadoop_bam.VCFFormat;
import java.util.ArrayList;
import java.util.List;
import static org.bgi.flexlab.gaea.tools.haplotypecaller.HaplotypeCallerEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY;

/**
 * Created by wangwl on 2017/5/18.
 */
public class HaplotypeCallerOptions extends GaeaOptions implements HadoopOptions {
    private final static String SOFTWARE_NAME = "HaplotypeCaller";
    private final static String SOFTWARE_VERSION = "0.1";

    
    /**
     * added in 2017.05.23 by Zhanghaorui
     */
    
    /**
     * Set of alleles to use in genotyping
     * 
     * private String alleles = none;
     */
    
    /**
     * dbSNP file
     * 
     * private String dbsnp = none;
     */
     
    /**
     * Output the active region to this IGV formatted file
     * 
     * private String activeRegionOut = NA;
     */
     
    /*
     * Output the raw activity profile results in IGV format
     * 
     * private String activityProofileOut = NA;
     */
     
    /**
     * Write debug assembly graph information to this file
     * 
     * private String graphOutput = NA;
     */
       
     
    /**
     * 	File to which variants should be written
     */
    private String out = "stdout";

    
    /**
     * Fraction of contamination to aggressively remove
     */
    private double contaminationFractionToFilter = 0.0;
    
    /**
     * Specifies how to determine the alternate alleles to use for genotyping
     */
    private String genotypingMode = "DISCOVERY";
    
    /**
     * One or more classes/groups of annotations to apply to variant calls
     */
    private String group = "[StandardAnnotation,StandardHCAnnotation]";
    
    /**
     * Heterozygosity value used to compute prior likelihoods for any locus
     */
    private double heterozygosity = 0.001;
    
    /**
     * Standard deviation of eterozygosity for SNP and indel calling.
     */
    private double heterozygosityStdev = 0.01;
    
    /**
     * Heterozygosity for indel calling
     */
    private double indelHeterozygosity = 1.25E-4;
    
    /**
     * Maximum reads in an active region
     */
    private int maxReadsInRegionPerSample = 10000;

    /**
     * Minimum base quality required to consider a base for calling
     */
    private byte minBaseQuality = 10;
    
    /**
     * Minimum number of reads sharing the same alignment start for each genomic location in an active region
     */
    private int minReadsPerAlignmentStart = 10;
    
    /**
     * Name of single sample to use from a multi-sample bam
     * 
     * private String sampleName = NA;
     */
    
    /**
     * Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).
     */
    private int samplePloidy = 2;
    
    /**
     * The minimum phred-scaled confidence threshold at which variants should be called
     */
    private double standardConfidenceForCalling = 10.0;
    
    /**
     * Annotate number of alleles observed
     */
    private boolean annotateNDA = false;
    
    /**
     * Use new AF model instead of the so-called exact model
     */
    private boolean useNewAFCalculator = false;
    
    /**
     * Use this interval list file as the active regions to process
     * 
     * private String activeRegionln = NA;
     */
    
    /**
     * 	Comparison VCF file
     * 
     * private String[] comp = "[]";
     */
    
    /**
     * File to which assembled haplotypes should be written
     * 
     * private String bamOutput = NA;
     */
    
    /**
     * Threshold for the probability of a profile state being active.
     */
    private double activeProbabilityThreshold = 0.002;
    
    /**
     * The active region extension; if not provided defaults to Walker annotated default
     * 
     * private String activeRegionExtension = NA;
     */
    
    /**
     * The active region maximum size; if not provided defaults to Walker annotated default
     * 
     * private String activeRegionMaxsize = NA;
     */
    
    /**
     * 	One or more specific annotations to apply to variant calls
     * 
     * private String[] annotation = "[]";
     */
    
    /**
     * Which haplotypes should be written to the BAM
     */
    private String bamWriterType = "CALLED_HAPLOTYPES";
    
    /**
     * The sigma of the band pass filter Gaussian kernel; if not provided defaults to Walker annotated default
     * 
     * private String bandPassSigma = NA;
     */
    
    /**
     * 	Contamination per sample
     * 
     * private String contaminationFractionPerSampleFile = NA;
     */
    
    /**
     * Mode for emitting reference confidence scores
     */
    private boolean emitRefConfidence = false;
    
    /**
     * One or more specific annotations to exclude
     * 
     * private String[] excludeAnnotation = "[]";
     */
    
    /**
     * Flat gap continuation penalty for use in the Pair HMM
     */
    private int gcpHMM = 10;
    
    /**
     * Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)
     */
    private int[] GVCFGQBands = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
    		                     23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
    		                     43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 70, 80, 90, 99};
   
    /**
     * The size of an indel to check for in the reference model
     */
    private int indelSizeToEliminateInRefModel = 10;
    
    /**
     * Input prior for calls
     * 
     * private String[] inputPrior = "[]";
     */
    
    /**
     * Kmer size to use in the read threading assembler
     */
    private int[]  kmerSize = {10,25};
    
    /**
     * Maximum number of alternate alleles to genotype
     */
    private int maxAlternateAlleles = 6;
    
    /**
     * Maximum number of genotypes to consider at any site
     */
    private int maxGenotypeCount = 1024;
    
    /**
     * Maximum number of PL values to output
     */
    private int maxNumPLValue = 100;
    
    /**
     * Maximum number of haplotypes to consider for your population
     */
    private int maxNumHaplotypesInPopulation = 128;
    
    /**
     * 	Maximum reads per sample given to traversal map() function
     */
    private int maxReadsInMemoryPerSample = 30000;
    
    /**
     * Maximum total reads given to traversal map() function
     */
    private int maxTotalReadsInMemory = 10000000;
    
    /**
     * 	Minimum length of a dangling branch to attempt recovery
     */
    private int minDanglingBranchLength = 4;
    
    /**
     * Minimum support to not prune paths in the graph
     */
    private int minPruning = 2;
    
    /**
     * Number of samples that must pass the minPruning threshold
     */
    private int numPruningSamples = 1;
    
    /**
     * Which type of calls we should output
     */
    private HaplotypeCallerEngine.OUTPUT_MODE outputMode = EMIT_VARIANTS_ONLY;
    
    /**
     * The PCR indel model to use
     */
    private String pcrIndleModel = "CONSERVATIVE";
    
    /**
     * The global assumed mismapping rate for reads
     */
    private int phredScaledGlobalReadMismappingRate = 45;
    
    /**
     * 	Allow graphs that have non-unique kmers in the reference
     */
    private boolean allowNonUniqueKmersInRef = false;
    
    /**
     * 	Annotate all sites with PLs
     */
    private boolean allSitePLs = false;
    
    /**
     * 1000G consensus mode
     */
    private boolean consensus = false;
    
    /**
     * Print out very verbose debug information about each triggering active region
     */
    private boolean debug = false;
    
    /**
     * Don't skip calculations in ActiveRegions with no variants
     */
    private boolean disableOptimizations = false;
    
    /**
     * Disable physical phasing
     */
    private boolean doNotRunPhysicalPhasing = false;
    
    /**
     * Disable iterating over kmer sizes when graph cycles are detected
     */
    private boolean dontIncreaseKmerSizesForCycles = false;
    
    /**
     * If specified, we will not trim down the active region from the full 
     * region (active + extension) to just the active interval for genotyping
     */
    private boolean dontTrimActiveRegions = false;
    
    /**
     * Do not analyze soft clipped bases in the reads
     */
    private boolean dontUseSoftClippedBases = false;
    
    /**
     * Emit reads that are dropped for filtering, trimming, realignment failure
     */
    private boolean emitDroppedReads = false;
    
    /**
     * If provided, all bases will be tagged as active
     */
    private boolean forceActive = false;
    
    /**
     * Use additional trigger on variants found in an external alleles file
     */
    private boolean useAllelesTrigger = false;	
    
    /**
     * Use the contamination-filtered read maps for the purposes of annotating variants
     */
    private boolean useFilteredReadsForAnnotations = false;
    
    /**
     * end here
     */
    
    
    public HaplotypeCallerOptions() {
        addOption("i", "input", true, "Input file containing sequence data (BAM or CRAM)");
        addOption("I", "is_sam_input", false, "the input is in SAM format.");
        addOption("o", "output", true, "directory to which variants should be written in HDFS written format.");
        addOption("O", "outputMode", true, "output mode :EMIT_VARIANTS_ONLY, EMIT_ALL_CONFIDENT_SITES, EMIT_ALL_SITES");
        addOption("B","is_bcf_output", false, "output variant in BCF format." );
        addOption("b", "bed_file", true, "only variant in this region will be called.");
        addOption("r", "reference", true, "Gaea Indexed reference list for memory sharing.");
        addOption("ARO", "activeRegionOut", true, "Output the active region to this IGV formatted file.");
        addOption("APO", "activityProfileOut", true, "Output the raw activity profile results in IGV format");
        addOption("graph", "graphOutput", true, "Write debug assembly graph information to this file");
        addOption("o","out",true,"File to which variants should be written");
        addOption("A", "annotation", true, "One or more specific annotations to apply to variant calls, the tag is separated by \',\'");
        addOption("contamination", "contamination_fraction_to_filter", true, "Fraction of contamination to aggressively remove.");
        addOption("gt_mode", "genotyping_mode", true, "Specifies how to determine the alternate alleles to use for genotyping");
        addOption("G", "group", true, "One or more classes/groups of annotations to apply to variant calls. The single value 'none' removes the default group, the group tag is separated by \',\'");
        addOption("hets", "heterozygosity", true, "Heterozygosity value used to compute prior likelihoods for any locus.");
        addOption("heterozygosityStandardDeviation", "heterozygosity_stdev", true, "\tStandard deviation of eterozygosity for SNP and indel calling.");
        addOption("indelHeterozygosity", "indel_heterozygosity", true, "Heterozygosity for indel calling.");
        addOption("maxReads", "maxReadsInRegionPerSample", true, "Maximum reads in an active region");
        addOption("mbq", "min_base_quality_score", true, "Minimum base quality required to consider a base for calling.");
        addOption("minReadsPerAlignStart", "minReadsPerAlignmentStart", true, "Minimum number of reads sharing the same alignment start for each genomic location in an active region");
        addOption("sn", "sample_name", true, "Name of single sample to use from a multi-sample bam");
        addOption("ploidy", "sample_ploidy", true, "Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).");
        addOption("standCallConf", "standardConfidenceForCalling", true, "standard confidence for calling");
        addOption("nda", "annotateNDA", false, "Annotate number of alleles observed");
        addOption("newQual", "useNewAFCalculator", false, "Use new AF model instead of the so-called exact model");
        addOption("AR", "activeRegionIn", true, "Use this interval list file as the active regions to process");
        addOption("comp", "comp", true, "Comparison VCF file");
        addOption("bamout", "bamOutput", true, "File to which assembled haplotypes should be written");
        addOption("ActProbThresh", "activeProbabilityThreshold", true, "Threshold for the probability of a profile state being active.");
        addOption("actRegionExtension", "activeRegionExtension", true, "The active region extension; if not provided defaults to Walker annotated default");
        addOption("actRegionMax", "activeRegionMaxSize", true, "The active region maximum size; if not provided defaults to Walker annotated default");
        addOption("A", "annotation", true, "One or more specific annotations to apply to variant calls");
        addOption("bamWriterType", "bamWriterType", true, "Which haplotypes should be written to the BAM");
        addOption("bandPassSigma", "bandPassSigma", true, "The sigma of the band pass filter Gaussian kernel; if not provided defaults to Walker annotated default");
        addOption("contaminationFile", "contamination_fraction_per_sample_file", true, "Contamination per sample");
        addOption("ERC", "emitRefConfidence", false, "Mode for emitting reference confidence scores");
        addOption("XA", "excludeAnnotation", true, "One or more specific annotations to exclude");
        addOption("gcpHMM", "gcpHMM", true, "Flat gap continuation penalty for use in the Pair HMM");
        addOption("GQB", "GVCFGQBands", true, "Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)");
        addOption("ERCIS", "indelSizeToEliminateInRefModel", true, "The size of an indel to check for in the reference model");
        addOption("inputPrior", "input_prior", true, "Input prior for calls");
        addOption("kmerSize", "kmerSize", true, "Kmer size to use in the read threading assembler");
        addOption("maxAltAlleles", "max_alternate_alleles", true, "Maximum number of alternate alleles to genotype");
        addOption("maxGT", "max_genotype_count", true, "Maximum number of genotypes to consider at any site");
        addOption("maxNumPLValues", "max_num_PL_values", true, "Maximum number of PL values to output");
        addOption("maxNumHaplotype", "maxNumHaplotypesInPopulation", true, "Maximum number of haplotypes to consider for your population");
        addOption("maxReadsPerSample", "maxReadsInMemoryPerSample", true, "Maximum reads per sample given to traversal map() function");
        addOption("maxTotalReads", "maxTotalReadsInMemory", true, "Maximum total reads given to traversal map() function");
        addOption("minDangling", "minDanglingBranchLength", true, "Minimum length of a dangling branch to attempt recovery");
        addOption("minPruning", "minPruning", true, "Minimum support to not prune paths in the graph");
        addOption("numPruningSamples", "numPruningSamples", true, "Number of samples that must pass the minPruning threshold");
        addOption("out_mode", "output_mode", true, "Which type of calls we should output");
        addOption("pcrModel", "pcr_indel_model", true, "The PCR indel model to use");
        addOption("globalMAPQ", "phredScaledGlobalReadMismappingRate", true, "The global assumed mismapping rate for reads");
        addOption("allowNonUniqueKmers", "allowNonUniqueKmersInRef", false, "Allow graphs that have non-unique kmers in the reference");
        addOption("allSitePLs", "allSitePLs", false, "Annotate all sites with PLs");
        addOption("consensus", "consensus", false, "1000G consensus mode");
        addOption("debug", "debug", false, "Print out very verbose debug information about each triggering active region");
        addOption("disableOptimizations", "disableOptimizations", false, "Don't skip calculations in ActiveRegions with no variants");
        addOption("doNotRunPhysicalPhasing", "doNotRunPhysicalPhasing", false, "Disable physical phasing");
        addOption("dontIncreaseKmerSizesForCycles", "dontIncreaseKmerSizesForCycles", false, "Disable iterating over kmer sizes when graph cycles are detected");
        addOption("dontTrimActiveRegions", "dontTrimActiveRegions", false, "If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping");
        addOption("dontUseSoftClippedBases", "dontUseSoftClippedBases", false, "Do not analyze soft clipped bases in the reads");
        addOption("edr", "emitDroppedReads", false, "Emit reads that are dropped for filtering, trimming, realignment failure");
        addOption("forceActive", "forceActive", false, "If provided, all bases will be tagged as active");
        addOption("allelesTrigger", "useAllelesTrigger", false, "Use additional trigger on variants found in an external alleles file");
        addOption("useFilteredReadsForAnnotations", "useFilteredReadsForAnnotations", false, "Use the contamination-filtered read maps for the purposes of annotating variants");
        addOption("S", "single_sample_mode", false, "will call genotype and variant for each sample separately");
        addOption("R", "reducer", true, "reducer numbers");
        addOption("W", "window_size", true, "window size that sharding the data.");
        addOption("h", "help", false, "print help information.");

        FormatHelpInfo(SOFTWARE_NAME, SOFTWARE_VERSION);
    }

    @Override
    public void parse(String[] args) {
        try {
            cmdLine = parser.parse(options, args);
        } catch (ParseException e) {
            printHelpInfotmation(SOFTWARE_NAME);
            throw new RuntimeException(e);
        }

        if(args.length == 0 || getOptionBooleanValue("h", false)) {
            printHelpInfotmation(SOFTWARE_NAME);
            System.exit(1);
        }
        input = getOptionValue("i", null);
        samFormat = getOptionBooleanValue("I", false);
        output = getOptionValue("o", null);
        bcfFormat = getOptionBooleanValue("B", false);
        bedRegionFile = getOptionValue("b", null);
        reference = getOptionValue("r", null);
        String annotationTags = getOptionValue("A", null);
        if(annotationTags != null) {
            for(String tag : annotationTags.split(","))
            annotations.add(tag);
        }
        annotationGroups.add("Standard");
        String annotationGroupTags = getOptionValue("G", null);
        if(annotationGroupTags != null) {
            if(annotationGroupTags.equals("none")) {
                annotationGroups.clear();
                annotationGroups.add("none");
            } else {
                for (String tag : annotationGroupTags.split(","))
                    annotationGroups.add(tag);
            }
        }
        contaminationFraction = getOptionDoubleValue("C1", DEFAULT_CONTAMINATION_FRACTION);
        heterozygosity = getOptionDoubleValue("hets", VariantCallingEngine.HUMAN_SNP_HETEROZYGOSITY);
        heterozygosityStandardDeviation = getOptionDoubleValue("heterozygosityStandardDeviation", 0.01);
        indelHeterozygosity = getOptionDoubleValue("indelHeterozygosity", 1.0/8000);
        maxDeletionFraction = getOptionDoubleValue("deletions", 0.05);
        noCapBaseQualsAtMappingQual = getOptionBooleanValue("C", false);
        minBaseQuality = getOptionByteValue("mbq", (byte)17);
        minMappingQuality = getOptionShortValue("mmq", (short)17);
        minIndelCountForGenotyping = getOptionIntValue("minIndelCnt", 5);
        minIndelFractionPerSample = getOptionDoubleValue("minIndelFrac", 0.25);
        samplePloidy = getOptionIntValue("ploidy", GaeaVariantContextUtils.DEFAULT_PLOIDY);
        standardConfidenceForCalling = getOptionDoubleValue("standCallConf", 10.0);
        standardConfidenceForEmitting = getOptionDoubleValue("standEmitConf", 30);
        indelGapContinuationPenalty = getOptionByteValue("indelGCP", (byte) 10);
        indelGapOpenPenalty = getOptionByteValue("indelGOP", (byte) 45);
        maxAlternateAlleles = getOptionIntValue("maxAltAlleles", 6);
        maxGenotypeCount = getOptionIntValue("maxGT", 1024);
        maxNumPLValues = getOptionIntValue("maxNumPLValues", 100);
        singleSampleMode = getOptionBooleanValue("S", false);
        reducerNumber = getOptionIntValue("R", 30);
        windowSize = getOptionIntValue("W", 100000);
        pcr_error = getOptionDoubleValue("pcrError", SNPGenotypeLikelihoodCalculator.DEFAULT_PCR_ERROR_RATE);

        try {
            gtlcalculators = GenotypeLikelihoodCalculator.Model.valueOf(getOptionValue("glm", "SNP"));
        } catch (Exception e) {
            throw new UserException.BadArgumentValueException("glm", e.getMessage());
        }
        try {
            outputMode = VariantCallingEngine.OUTPUT_MODE.valueOf(getOptionValue("O", "EMIT_VARIANTS_ONLY"));
        } catch (Exception e) {
            throw new UserException.BadArgumentValueException("O", e.getMessage());
        }
        try {
            AFmodel = AFCalcFactory.Calculation.valueOf(getOptionValue("AFmodel", "EXACT_INDEPENDENT"));
        } catch (Exception e) {
            throw new UserException.BadArgumentValueException("AFmodel", e.getMessage());
        }
        try {
            pairHmmImplementation = PairHMM.HMM_IMPLEMENTATION.valueOf(getOptionValue("pairHMM", "ORIGINAL"));
        } catch (Exception e) {
            throw new UserException.BadArgumentValueException("pairHMM", e.getMessage());
        }

        check();

        minDepth = getOptionIntValue("D", 4);
        annotateNumberOfAllelesDiscovered = getOptionBooleanValue("numAlleleDis", false);
        noCapBaseQualsAtMappingQual = getOptionBooleanValue("C", false);

    }

    private void check() {
        if(input == null)
            throw new UserException.BadArgumentValueException("i", "input directory or file is not assigned.");

        if(output == null)
            throw new UserException.BadArgumentValueException("o", "output directory is not assigned.");

        if(reference == null)
            throw new UserException.BadArgumentValueException("r", "reference can not be null.");

        if(reducerNumber <= 0 ) {
            throw new UserException.BadArgumentValueException("R", "reducer number can not be less than 1.");
        }
    }

    @Override
    public void setHadoopConf(String[] args, Configuration conf) {
        conf.setStrings("args", args);
    }

    @Override
    public void getOptionsFromHadoopConf(Configuration conf) {
        String[] args = conf.getStrings("args");
        this.parse(args);
    }

    public GenotypeLikelihoodCalculator.Model getGtlcalculators() {
        return gtlcalculators;
    }

    public void setGtlcalculators(GenotypeLikelihoodCalculator.Model gtlcalculators) {
        this.gtlcalculators = gtlcalculators;
    }

    public SAMFormat getInputFormat(){
        if(samFormat)
            return SAMFormat.SAM;
        return SAMFormat.BAM;
    }

    public VCFFormat getOuptputFormat() {
        if(bcfFormat)
            return VCFFormat.BCF;
        return VCFFormat.VCF;
    }

    public int getMinDepth() {
        return minDepth;
    }

    public void setMinDepth(int minDepth) {
        this.minDepth = minDepth;
    }

    public boolean isCapBaseQualsAtMappingQual() {
        return !noCapBaseQualsAtMappingQual;
    }

    public byte getMinBaseQuality() {
        return minBaseQuality;
    }

    public short getMinMappingQuality() {
        return minMappingQuality;
    }

    public VariantCallingEngine.OUTPUT_MODE getOutputMode() {
        return outputMode;
    }

    public int getMinIndelCountForGenotyping() {
        return minIndelCountForGenotyping;
    }

    public double getMinIndelFractionPerSample() {
        return minIndelFractionPerSample;
    }

    public double getContaminationFraction() {
        return contaminationFraction;
    }

    public int getMaxAlternateAlleles() {
        return maxAlternateAlleles;
    }

    public int getSamplePloidy() {
        return samplePloidy;
    }

    public AFCalcFactory.Calculation getAFmodel() {
        return AFmodel;
    }

    public void setAFmodel(AFCalcFactory.Calculation AFmodel) {
        this.AFmodel = AFmodel;
    }

    public double getStandardConfidenceForCalling() {
        return standardConfidenceForCalling;
    }

    public double getStandardConfidenceForEmitting() {
        return standardConfidenceForEmitting;
    }

    public double getHeterozygosity() {
        return heterozygosity;
    }

    public double getIndelHeterozygosity() {
        return indelHeterozygosity;
    }

    public boolean isAnnotateNumberOfAllelesDiscovered() {
        return annotateNumberOfAllelesDiscovered;
    }

    public boolean isSingleSampleMode() {
        return singleSampleMode;
    }

    public String getInput() {
        return input;
    }

    public String getOutput() {
        return output;
    }

    public String getBAMHeaderOutput(){
        return output+"bamHeader";
    }

    public String getVCFHeaderOutput() {
        return output + "vcfHeader";
    }

    public String getReference() {
        return  reference;
    }

    public int getReducerNumber() {
        return reducerNumber;
    }

    public int getWindowSize() {
        return windowSize;
    }

    public List<String> getAnnotations() {
        return annotations;
    }

    public List<String> getAnnotationGroups() {
        return annotationGroups;
    }

    public String getBedRegionFile() {
        return bedRegionFile;
    }

    public double getPcr_error() {
        return pcr_error;
    }

    public byte getIndelGapOpenPenalty() {
        return indelGapOpenPenalty;
    }

    public byte getIndelGapContinuationPenalty() {
        return indelGapContinuationPenalty;
    }

    public PairHMM.HMM_IMPLEMENTATION getPairHmmImplementation() {
        return pairHmmImplementation;
    }
}

