package org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer.Context;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.structure.positioninformation.depth.PositionDepth;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.region.RegionChromosomeInformation;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion;
import org.bgi.flexlab.gaea.data.structure.region.TargetRegion;
import org.bgi.flexlab.gaea.data.structure.region.SingleRegion.Regiondata;
import org.bgi.flexlab.gaea.data.structure.region.statistic.BedSingleRegionStatistic;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.RegionReport.Sex;
import org.bgi.flexlab.gaea.tools.mapreduce.bamqualitycontrol.BamQualityControlOptions;
import org.bgi.flexlab.gaea.util.SamRecordDatum;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;

public class RegionResultReport extends ResultReport{
	
	private TargetRegion region;
	
	private SingleRegion cnvRegion;
	
	private SingleRegion bedRegion;
	
	private SingleRegion genderRegion;
	
	private RegionReport regionReport;
	
	private BedSingleRegionReport bedSingleRegionReport;
	
	private BedSingleRegionReport genderSingleRegionReport;
	
	private CNVDepthReport cnvDepthReport;
	
	private Map<String, Integer> sampleLaneSize = new HashMap<>();
	
	private Configuration conf;
					
	public RegionResultReport(BamQualityControlOptions options, Configuration conf) throws IOException {
		super(options);
		region = new TargetRegion();
		cnvRegion = new SingleRegion();
		if (options.getRegion() != null) {
			region.parseRegion(options.getRegion(), true);
		}
		if (options.getBedfile() != null) {
			region.parseBedFileFromHDFS(options.getBedfile(), true);
			bedRegion = new SingleRegion();
			bedRegion.parseRegionsFileFromHDFS(options.getBedfile(), true, 0);
			genderRegion = new SingleRegion();
			genderRegion.parseRegionsFileFromHDFS(options.getBedfile(), true, 0);
		}
		this.conf = conf;
		if(options.isCnvDepth())
			initSampleLaneSize(conf);
	}
	
	public void initReports(String sample) throws IOException {
		super.initReports();
		regionReport = new RegionReport(region);
		if(options.getBedfile() != null){
			bedSingleRegionReport = new BedSingleRegionReport(bedRegion);
			genderSingleRegionReport = new BedSingleRegionReport(genderRegion);
		}
			
		if(options.isCnvDepth() && options.getCnvRegion() != null) {
			cnvRegion = new SingleRegion();
			cnvRegion.parseRegionsFileFromHDFS(options.getCnvRegion(), true, 0);
			cnvDepthReport = new CNVDepthReport(sampleLaneSize.get(sample),cnvRegion);
		}
	}
	
	private void initSampleLaneSize(Configuration conf) {
		SAMFileHeader mFileHeader = SamHdfsFileHeader.getHeader(conf);
		for(SAMReadGroupRecord rg : mFileHeader.getReadGroups()) {
			if(sampleLaneSize.containsKey(rg.getSample())) {
				int tmp = sampleLaneSize.get(rg.getSample());
				sampleLaneSize.put(rg.getSample(), tmp + 1);
			} else {
				sampleLaneSize.put(rg.getSample(), 1);
			}
		}
	}
	
	@Override
	public boolean mappedReport(SamRecordDatum datum, String chrName, Context context) {
		super.mappedReport(datum, chrName, context);
		regionReport.constructMapReport(chrName, datum);
		return true;
	}
	
	@Override
	public void constructDepthReport(PositionDepth pd, int i, String chrName, long pos) {
		int depth = pd.getPosDepth(i);
		int noPCRdepth = pd.getRMDupPosDepth(i);
		if(region.isPositionInRegion(chrName, pos)) {
			super.regionCoverReport(depth, noPCRdepth);
			if(options.isOutputUnmapped() && depth != 0 )
				unmappedReport.updateUnmappedSites(pos, unmappedReport.getUnmappedSites(chrName));
		}
		if(options.isCnvDepth() && cnvRegion.posInRegion(chrName, (int) (pos)) != -1) {
			cnvDepthReport.add(chrName, (int)pos, pd.getLaneDepth(i));
		}
		
		regionReport.constructDepthReport(pd, i, chrName, pos);
	}
	
	@Override 
	public void singleRegionReports(String chrName, long winStart, int winSize , PositionDepth pd) {
		super.singleRegionReports(chrName, winStart, winSize, pd);
		if(options.getBedfile() != null) {
			bedSingleRegionReport.getStatisticString(chrName, (int)winStart, winSize, pd.getNormalPosDepth(), "bed");
			genderSingleRegionReport.getStatisticString(chrName,(int) winStart, winSize, pd.getGenderPosDepth(), "gender");
		}
	}
	
	@Override
	public String toReducerString(String sample, String chrName, boolean unmappedRegion) {
		StringBuffer info = new StringBuffer();
		
		if(chrName == "-1") {
			info.append("sample:");
			info.append(sample);
			info.append("\n");
			info.append(basicReport.toReducerString());
			return info.toString();
		}
		info.append("sample:");
		info.append(sample);
		info.append("\n");
		info.append("chrName:");
		info.append(chrName);
		info.append("\n");
		info.append(basicReport.toReducerString());
		if(!unmappedRegion) {
			info.append(regionCoverReport.toReducerString());
			info.append(regionReport.toReducerString());
			if(options.isCnvDepth()) {
				//System.err.println("do cnv depth");
				info.append(cnvDepthReport.toReducerString());
			}
			info.append("insert size information:\n");
			insertSizeReportReducerString(info, insertSize);
			info.append("insert size information\n");
			
			info.append("insert size without dup information:\n");
			insertSizeReportReducerString(info, insertSizeWithoutDup);
			info.append("insert size without dup information\n");
			
			info.append("unmapped site information:\n");
			unmappedReport.toReducerString();
			info.append("unmapped site information\n");
			if(cnvSingleRegionReport != null)
				info.append(cnvSingleRegionReport.toReducerString());
			if(bedSingleRegionReport != null) 
				info.append(bedSingleRegionReport.toReducerString());
			if(genderSingleRegionReport != null) 
				info.append(genderSingleRegionReport.toReducerString());
		}
		return info.toString();
	}
	
	@Override
	public void parseReport(LineReader lineReader, Text line, ReferenceShare genome) throws IOException {
		super.parseReport(lineReader, line,genome);
		String lineString = line.toString();
		if(lineString.contains("Target Information")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				regionReport.parse(line.toString());
			}
		}
		if(lineString.startsWith("bed single Region Statistic")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				bedSingleRegionReport.parseReducerOutput(line.toString(), false);
			}
		}
		if(lineString.startsWith("bed part single Region Statistic")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				bedSingleRegionReport.parseReducerOutput(line.toString(), true);
			}
		}
		if(lineString.startsWith("gender single Region Statistic")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				bedSingleRegionReport.parseReducerOutput(line.toString(), false);
			}
		}
		if(lineString.startsWith("gender part single Region Statistic")) {
			if(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				bedSingleRegionReport.parseReducerOutput(line.toString(), true);
			}
		}
		if(lineString.startsWith("CNV Depth")) {
			while(lineReader.readLine(line) > 0 && line.getLength() != 0) {
				if(line.toString().contains("CNV Depth")) {
					break;
				}
				cnvDepthReport.add(line.toString());
			}
		}
	}
	
	@Override
	public void write(FileSystem fs, String sampleName) throws IOException {
		super.write(fs, sampleName);
		if(bedSingleRegionReport != null) {
			//chromosome info
			Map<Regiondata, BedSingleRegionStatistic> result = bedSingleRegionReport.getResult();
			Map<String, RegionChromosomeInformation> chrsInfo = new HashMap<String, RegionChromosomeInformation>();
			
			for(Regiondata regionData : result.keySet()) {
				if(chrsInfo.containsKey(regionData.getChrName())) {
					chrsInfo.get(regionData.getChrName()).add(regionData, result.get(regionData));
				} else {
					RegionChromosomeInformation regionChrInfo = new RegionChromosomeInformation();
					regionChrInfo.add(regionData, result.get(regionData));
					chrsInfo.put(regionData.getChrName(), regionChrInfo);
				}
			}
			StringBuffer regionChromosomeFilePath = new StringBuffer();
			regionChromosomeFilePath.append(options.getOutputPath());
			regionChromosomeFilePath.append("/");
			regionChromosomeFilePath.append(sampleName);
			regionChromosomeFilePath.append(".chromosome.txt");
			Path singleRegionPath = new Path(regionChromosomeFilePath.toString());
			FSDataOutputStream regionChromosomeWriter = fs.create(singleRegionPath);
			for(String chrName : chrsInfo.keySet()) {
				regionChromosomeWriter.write(chrsInfo.get(chrName).toString(chrName).getBytes());
			}
			regionChromosomeWriter.close();
		}
		
		if(genderSingleRegionReport != null) {
			Map<Regiondata, BedSingleRegionStatistic> result = genderSingleRegionReport.getResult();
			ArrayList<Double> depth = new ArrayList<Double>();
			
			for(Regiondata regionData : result.keySet()) {
				depth.add(result.get(regionData).getAverageDepth());
				ChromosomeInformationShare chrInfo = genome.getChromosomeInfo(regionData.getChrName());
				result.get(regionData).calRefGCrate(chrInfo, regionData, 0);
			}
			
			if(!options.isGenderPredict() || result.size() < 10) {
				System.err.println("[Warnning] input bed file have too little regions, abort gender cal!");
				regionReport.setPredictedGender(Sex.unKown);
			} else {
				Collections.sort(depth);
				double[] filterValue = new double[2];
				filterValue[0] = depth.get((int) (depth.size() * 0.02));
				filterValue[1] = depth.get((int) (depth.size() * 0.98));
//				GenderPredict gender = new GenderPredict(result, genderSingleRegionReport.getRegion(), filterValue);
//				regionReport.setPredictedGender(gender.predict());
			}
		}
		
		StringBuffer reportFilePath = new StringBuffer();
		StringBuffer info = new StringBuffer();
		reportFilePath.append(options.getOutputPath());
		reportFilePath.append("/");
		reportFilePath.append(sampleName);
		reportFilePath.append(".bam.report.txt");
		Path reportPath = new Path(reportFilePath.toString());
		FSDataOutputStream reportwriter = fs.create(reportPath);
		info.append(basicReport.toString());
		info.append(regionReport.toString(basicReport));
		reportwriter.write(info.toString().getBytes());
		reportwriter.close();
		
		if(cnvDepthReport != null)
			cnvDepthReport.toReport(options, fs, conf, sampleName);
	}
	
	public int getSampleLaneSize(String sample) {
		return sampleLaneSize.get(sample);
	}
	
	public BedSingleRegionReport getBedSingleRegionReport() {
		return bedSingleRegionReport;
	}
	
	public BedSingleRegionReport getGenderSingleRegionReport() {
		return genderSingleRegionReport;
	}
}
