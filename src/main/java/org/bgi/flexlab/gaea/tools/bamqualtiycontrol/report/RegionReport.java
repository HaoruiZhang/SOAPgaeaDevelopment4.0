package org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map.Entry;

import org.bgi.flexlab.gaea.data.structure.positioninformation.depth.PositionDepth;
import org.bgi.flexlab.gaea.data.structure.region.TargetRegion;
import org.bgi.flexlab.gaea.data.structure.region.FlankRegion;
import org.bgi.flexlab.gaea.data.structure.region.Region;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.SamRecordDatum;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.CounterProperty.BaseType;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.CounterProperty.Depth;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.CounterProperty.DepthType;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.CounterProperty.Interval;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.CounterProperty.ReadType;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.Tracker.BaseTracker;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.report.Tracker.ReadsTracker;

public class RegionReport{
	
	private TargetRegion region;
			
	private Sex gender;
		
	private BaseTracker bTracker;
	
	private ReadsTracker rTracker;
	
	public RegionReport(Region region) {
		bTracker = new Tracker.BaseTracker();
		rTracker = new Tracker.ReadsTracker();
		register();
		this.region = (TargetRegion) region;
	}
	
	public enum Sex{
		M,
		F,
		unKown
	}
	
	public void constructMapReport(String chrName, SamRecordDatum datum) {
		if (region.isReadInRegion(chrName, datum.getPosition(), datum.getEnd())) {
			rTracker.setTrackerAttribute(ReadType.MAPPED, Interval.TARGET);
			rTracker.setTrackerAttribute(ReadType.MAPPED, Interval.FLANK);
			if (datum.isUniqueAlignment()) {
				rTracker.setTrackerAttribute(ReadType.UNIQUE, Interval.TARGET);
			}
		} else {
			if (region.isReadInFlank(chrName, datum.getPosition(), datum.getEnd())) {
				rTracker.setTrackerAttribute(ReadType.MAPPED, Interval.FLANK);
			}
		}	
	}
	
	public void constructDepthReport(PositionDepth pd, int i, String chrName, long pos ) {
		int depth = pd.getPosDepth(i);
		if(depth != 0) {
			int noPCRdepth = pd.getRMDupPosDepth(i);
			if(region.isPositionInRegion(chrName, pos)) {
				bTracker.setTrackerAttribute(Interval.TARGET, DepthType.NORMAL.setDepth(depth), DepthType.WITHOUT_PCR.setDepth(noPCRdepth));
				if(chrName.equals("chrX") || chrName.equals("chrx")) {
					bTracker.setTrackerAttribute(Interval.CHRX, DepthType.NORMAL.setDepth(depth), DepthType.WITHOUT_PCR);
				}
				if(chrName.equals("chrY") || chrName.equals("chry")) {
					bTracker.setTrackerAttribute(Interval.CHRY, DepthType.NORMAL.setDepth(depth), DepthType.WITHOUT_PCR);
				}
			} 
			if(region.isPositionInFlank(chrName, pos)) 
				bTracker.setTrackerAttribute(Interval.FLANK, DepthType.NORMAL.setDepth(depth), DepthType.WITHOUT_PCR.setDepth(noPCRdepth));
		} else if(depth == 0) {
			if(region.isPositionInRegion(chrName, pos) && pd.isDeletionBaseWithNoConver(i)) {
				bTracker.setTrackerAttribute(Interval.TARGET, DepthType.NORMAL, DepthType.WITHOUT_PCR);
			}
		}
	}
	
	public String toReducerString() {
		StringBuffer regionString = new StringBuffer();
		regionString.append("Target Information:\n");
		for(Entry<String, ReadsCounter> counter : rTracker.getCounterMap().entrySet()) {
			regionString.append(counter.getKey());
			regionString.append(" ");
			regionString.append(counter.getValue().getReadsCount());
			regionString.append("\t");
		}
		for(Entry<String, BaseCounter> counter : bTracker.getCounterMap().entrySet()) {
			regionString.append(counter.getKey());
			regionString.append(" ");
			regionString.append(counter.getValue().getProperty());
			regionString.append("\t");
		}
		regionString.append("\n");
		return regionString.toString();
	}
	
	public String toString(BasicReport basicReport) {
		DecimalFormat df = new DecimalFormat("0.000");
		df.setRoundingMode(RoundingMode.HALF_UP);
		
		StringBuffer regionString = new StringBuffer();
		regionString.append("Target Information:\n");
		regionString.append("Mapped Reads Number in Target:\t");
		regionString.append(rTracker.getReadsCount(ReadType.MAPPED, Interval.TARGET));
		regionString.append("\nUniq Mapped Reads Number in Target:\t");
		regionString.append(rTracker.getReadsCount(ReadType.UNIQUE, Interval.TARGET));
		regionString.append("\nCapture specificity:\t");
		regionString.append(df.format(getCaptureSpecificity(basicReport)));
		regionString.append("%\nCapture Effiency:\t");
		regionString.append(df.format(getCaptureEffiency(basicReport)));
		regionString.append("%\nMapped Reads Number in Flanking:\t");
		regionString.append(rTracker.getReadsCount(ReadType.MAPPED, Interval.FLANK));
		regionString.append("\nTarget size:\t");
		regionString.append(region.getRegionSize());
		formateCoverageInfo(regionString, "Target", region, Interval.TARGET, DepthType.NORMAL);
		regionString.append("\nFlanking size:\t");
		FlankRegion fRegion = region.getFlankRegion();
		regionString.append(fRegion.getRegionSize());
		formateCoverageInfo(regionString, "Flanking", fRegion, Interval.FLANK, DepthType.NORMAL);
		formateCoverageInfo(regionString, "Target", region, Interval.TARGET, DepthType.WITHOUT_PCR);
		formateCoverageInfo(regionString, "Flanking", fRegion, Interval.FLANK, DepthType.WITHOUT_PCR);
		regionString.append("\nchrX region coverage:\t");
		regionString.append(df.format(100 * (bTracker.getBaseCount(Interval.CHRX, Depth.ABOVE_ZREO, DepthType.NORMAL)/(double)region.getChrSize("chrX"))));
		regionString.append("%\nchrX region depth:\t");
		regionString.append(df.format(getMeanDepth(Interval.CHRX)));
		regionString.append("\nchrY region coverage:\t");
		regionString.append(df.format(100 * (bTracker.getBaseCount(Interval.CHRY, Depth.ABOVE_ZREO, DepthType.NORMAL)/(double)region.getChrSize("chrY"))));
		regionString.append("%\nchrY region depth:\t");
		regionString.append(df.format(getMeanDepth(Interval.CHRY)));
		regionString.append("\npredicted gender:\t");
		regionString.append(gender.toString());
		regionString.append("\n");
		return regionString.toString();
	}

	private void formateCoverageInfo(StringBuffer regionString, String location, Region region, Interval interval, DepthType type) {
		DecimalFormat df = new DecimalFormat("0.000");
		df.setRoundingMode(RoundingMode.HALF_UP);
		String rmdupSuffix = type == DepthType.WITHOUT_PCR ? "[RM DUP]" : "";
		regionString.append(String.format("\n%s coverage%s:\t", location, rmdupSuffix));
		regionString.append(df.format(getCoverage(interval, Depth.ABOVE_ZREO, type, region.getRegionSize())));
		regionString.append(String.format("%\n%s coverage > 4X percentage%s:\t", location, rmdupSuffix));
		regionString.append(df.format(100 * (bTracker.getBaseCount(interval, Depth.FOURX, type) /(double)region.getRegionSize())));
		regionString.append(String.format("%\n%s coverage > 10X percentage%s:\t", location, rmdupSuffix));
		regionString.append(df.format(100 * (bTracker.getBaseCount(interval, Depth.TENX, type)/(double)region.getRegionSize())));
		regionString.append(String.format("%\n%s coverage > 20X percentage%s:\t", location, rmdupSuffix));
		regionString.append(df.format(100 * (bTracker.getBaseCount(interval, Depth.TWENTYX, type)/(double)region.getRegionSize())));
		regionString.append(String.format("%\n%s coverage > 30X percentage%s:\t", location, rmdupSuffix));
		regionString.append(df.format(100 * (bTracker.getBaseCount(interval, Depth.THIRTYX, type)/(double)region.getRegionSize())));
		regionString.append(String.format("%\n%s coverage > 50X percentage%s:\t", location, rmdupSuffix));
		regionString.append(df.format(100 * (bTracker.getBaseCount(interval, Depth.FIFTYX, type)/(double)region.getRegionSize())));
		regionString.append(String.format("%\n%s coverage > 100X percentage%s:\t", location, rmdupSuffix));
		regionString.append(df.format(100 * (bTracker.getBaseCount(interval, Depth.HUNDREDX, type)/(double)region.getRegionSize())));
		regionString.append(String.format("%\n%s Mean Depth:\t", location));
		regionString.append(df.format(getMeanDepth(interval)));
	}
	
	public void register() {
		bTracker.register(createBaseCounters());
		rTracker.register(createReadsCounters());
	}
	
	public BaseTracker getBaseTracker() {
		return bTracker;
	}
	
	public ReadsTracker getReadsTracker() {
		return rTracker;
	}
	
	public List<BaseCounter> createBaseCounters() {
		List<BaseCounter> observers = new ArrayList<>();
		for(Interval region : Interval.values()) {
			switch (region) {
			case CHRX:
			case CHRY:
				observers.add(new BaseCounter(region, Depth.TOTALDEPTH, DepthType.NORMAL));
				observers.add(new BaseCounter(region, Depth.ABOVE_ZREO, DepthType.NORMAL));
			default:
				for(Depth depth : Depth.values())
					for(DepthType type : DepthType.values()) {
						observers.add(new BaseCounter(region, depth, type));
						};
			}
		}
		return observers;
	}
	
	public List<ReadsCounter> createReadsCounters() {
		List<ReadsCounter> counters = new ArrayList<>();
		Collections.addAll(counters, new ReadsCounter(ReadType.MAPPED, Interval.TARGET),
									 new ReadsCounter(ReadType.MAPPED, Interval.FLANK),
									 new ReadsCounter(ReadType.UNIQUE, Interval.TARGET));
		return counters;
	}
	
	public long getTotalDepth(Interval region, Depth depth, DepthType type) {
		return bTracker.getTotalDepth(region,depth,type);
	}
	
	public double getMeanDepth(Interval region) {
		long regionDepth = bTracker.getTotalDepth(region, Depth.TOTALDEPTH, DepthType.NORMAL);
		long baseCount = bTracker.getBaseCount(region, Depth.ABOVE_ZREO, DepthType.NORMAL);
		return regionDepth * 1.0 / baseCount;
	}
	
	public void setPredictedGender(Sex gender) {
		this.gender = gender;
	}
	
	public Sex getPredictedGender() {
		return this.gender;
	}
	
	public double getCoverage(Interval region, Depth depth, DepthType type, int size) {
		return (100 * (bTracker.getBaseCount(region, depth, type)/(double)size));
	}
	
	public double getCaptureSpecificity(BasicReport basicReport) {
		return 100 * (rTracker.getReadsCount(ReadType.UNIQUE, Interval.TARGET) / (double) basicReport.getReadsTracker().getReadsCount(ReadType.UNIQUE));
	}
	
	public double getCaptureEffiency(BasicReport basicReport) {
		return 100 * (bTracker.getTotalDepth(Interval.TARGET, Depth.TOTALDEPTH, DepthType.NORMAL) / (double) basicReport.getBaseTracker().getBaseCount(BaseType.TOTALBASE));
	}
	
}
